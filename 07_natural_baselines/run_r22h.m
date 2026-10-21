function result=run_r22h(root,stage)
    pkg=fileparts(mfilename('fullpath'));pkg=fileparts(pkg);
    audit=fullfile(root,'pufb_otfs_v2','audit','r2_2h_natural_baselines');
    output=fullfile(root,'pufb_otfs_v2','results','r2_2h_natural_baselines');
    r22d_audit=fullfile(pkg,'03_axis_specific_pufb');
    r22d_results=fullfile(root,'pufb_otfs_v2','results','r2_2d_axis_specific');
    r22g_audit=fullfile(root,'pufb_otfs_v2','audit','r2_2g_factorial_ablation');
    r22g_results=fullfile(root,'pufb_otfs_v2','results','r2_2g_factorial_ablation');
    r22a_audit=fullfile(pkg,'02_shared_lib');
    r20a_audit=fullfile(pkg,'01_channel_generator');
    r10r_audit=fullfile(root,'pufb_otfs_v2','audit','r1_0r_wotfs_reproduction');
    addpath(audit,r22d_audit,r22g_audit,r22a_audit,r20a_audit,r10r_audit);
    p=jsondecode(fileread(fullfile(audit,'r22h_preregistration.json')));
    switch upper(stage)
        case 'PREPARE'
            result=prepare_r22h(root,audit,output,r22d_results,r22g_results,p);
        case 'TEST'
            result=evaluate_r22h(root,audit,output,r22d_results,r22g_results,p);
        otherwise
            error('R22H:Stage','Use PREPARE or TEST.');
    end
end

function result=prepare_r22h(root,audit,output,r22d_results,r22g_results,p)
    if isfolder(output)&&isfile(fullfile(output,'r22h_execution_status.json'))
        result=refresh_pretest_freeze(root,audit,output,r22d_results,r22g_results,p);return
    end
    partial_resume=isfolder(output);
    if partial_resume
        listing=dir(output);existing=string({listing(~[listing.isdir]).name});
        allowed=["r22h_baseline_definitions.md";"r22h_preregistration.json"; ...
            "r22h_preregistration_amendment_01.json";"r22h_preregistration_amendment_02.json"; ...
            "r22h_support_hashes.csv"];
        assert(all(ismember(existing,allowed))&&~isfile(fullfile(output,'r22h_execution_status.json')), ...
            'R22H:OutputExists','Refusing to overwrite an existing or non-partial R2.2H results directory.');
        assert(strcmp(r22g_math('sha_file',fullfile(output,'r22h_preregistration.json')), ...
            r22g_math('sha_file',fullfile(audit,'r22h_preregistration.json'))));
        assert(strcmp(r22g_math('sha_file',fullfile(output,'r22h_baseline_definitions.md')), ...
            r22g_math('sha_file',fullfile(audit,'r22h_baseline_definitions.md'))));
    end
    hashes=verify_frozen_inputs(root,r22d_results,r22g_results,p);
    definition_check=check_baseline_definitions();
    train_manifest=readtable(fullfile(r22d_results,'r22d_train_manifest.csv'),'TextType','string');
    val_manifest=readtable(fullfile(r22d_results,'r22d_val_manifest.csv'),'TextType','string');
    test_manifest=readtable(fullfile(r22d_results,'r22d_test_manifest.csv'),'TextType','string');
    validate_manifest(train_manifest,200,121100000,40);
    validate_manifest(val_manifest,100,121200000,20);
    validate_manifest(test_manifest,1000,121300000,200);
    [names,definitions,U,alpha_nu,alpha_tau,frozen_masks]=method_specs(r22g_results);
    train=load_partition(r22d_results,train_manifest,'train');
    masks=cell(1,5);retained=zeros(5,1);support_hash=strings(5,1);
    for method=1:5
        energy=zeros(288);
        for channel=1:size(train,3)
            H=U{method}'*train(:,:,channel)*U{method};energy=energy+abs(H).^2;
        end
        energy=energy/size(train,3);[~,order]=sort(energy(:),'descend');mask=false(288);mask(order(1:p.support_protocol.K))=true;
        masks{method}=mask;retained(method)=sum(energy(mask))/sum(energy,'all');
        support_hash(method)=string(r22g_math('support_hash',mask));
    end
    supports=r22g_math('read',fullfile(r22g_results,'r22g_support_hashes.csv'));
    ir=supports(supports.Method=="I-R",:);hh=supports(supports.Method=="H-H",:);
    assert(height(ir)==1&&height(hh)==1);
    assert(isequal(masks{5},frozen_masks{5})&&support_hash(5)==string(ir.MaskSHA256), ...
        'R22H:IndependentSupport','Frozen R2.2G I-R transform/support did not reproduce.');
    assert(isequal(masks{4},frozen_masks{4})&&support_hash(4)==string(hh.MaskSHA256), ...
        'R22H:HaarSupport','R2.2G Haar support did not reproduce from R2.2D TRAIN.');
    if ~partial_resume,mkdir(output);end
    source_files={'r22h_preregistration.json','r22h_baseline_definitions.md', ...
        'r22h_preregistration_amendment_01.json','r22h_preregistration_amendment_02.json', ...
        'r22h_preregistration_amendment_03.json'};
    for source_index=1:numel(source_files)
        destination=fullfile(output,source_files{source_index});
        if ~isfile(destination),copyfile(fullfile(audit,source_files{source_index}),destination);end
    end
    support_table=table(names(:),definitions(:),alpha_nu(:),alpha_tau(:),repmat(p.support_protocol.K,5,1), ...
        support_hash,retained,'VariableNames',{'Representation','Definition','AlphaNu','AlphaTau','K','MaskSHA256','TrainRetainedEnergyFraction'});
    support_path=fullfile(output,'r22h_support_hashes.csv');
    if isfile(support_path)
        prior_support=r22g_math('read',support_path);
        numeric_columns={'AlphaNu','AlphaTau','K','TrainRetainedEnergyFraction'};
        numeric_match=true;
        for column=1:numel(numeric_columns)
            prior_values=prior_support.(numeric_columns{column});current_values=support_table.(numeric_columns{column});
            numeric_match=numeric_match&&all(abs(prior_values-current_values)<1e-14| ...
                (isnan(prior_values)&isnan(current_values)));
        end
        text_match=isequal(prior_support.Representation,support_table.Representation) ...
            &&isequal(prior_support.Definition,support_table.Definition) ...
            &&isequal(prior_support.MaskSHA256,support_table.MaskSHA256);
        assert(height(prior_support)==height(support_table)&&text_match&&numeric_match, ...
            'R22H:PartialSupportMismatch', ...
            'Existing partial TRAIN supports differ; refusing to overwrite.');
    else
        writetable(support_table,support_path);
    end
    overlap_rows=cell(5,6);independent=masks{5};
    for method=1:5
        intersection=nnz(independent&masks{method});union_count=nnz(independent|masks{method});
        overlap_rows(method,:)={names(method),intersection,intersection/p.support_protocol.K, ...
            intersection/union_count,support_hash(method),support_hash(5)};
    end
    writetable(cell2table(overlap_rows,'VariableNames',{'Representation','IntersectionCount', ...
        'IntersectionOverK','Jaccard','BaselineMaskSHA256','IndependentMaskSHA256'}), ...
        fullfile(output,'r22h_support_overlap.csv'));
    save(fullfile(output,'r22h_frozen_transforms.mat'),'U','masks','names','definitions','alpha_nu','alpha_tau','-v7.3');
    input_hashes=table(hashes.Paths,hashes.SHA256,'VariableNames',{'Path','SHA256'});
    writetable(input_hashes,fullfile(output,'r22h_input_hashes.csv'));
    implementation=make_implementation_hashes(root,audit);
    writetable(implementation,fullfile(output,'r22h_implementation_hashes.csv'));
    prereg_hash=string(r22g_math('sha_file',fullfile(audit,'r22h_preregistration.json')));
    definitions_hash=string(r22g_math('sha_file',fullfile(audit,'r22h_baseline_definitions.md')));
    amendment_hash=string(r22g_math('sha_file',fullfile(audit,'r22h_preregistration_amendment_01.json')));
    amendment_02_hash=string(r22g_math('sha_file',fullfile(audit,'r22h_preregistration_amendment_02.json')));
    amendment_03_hash=string(r22g_math('sha_file',fullfile(audit,'r22h_preregistration_amendment_03.json')));
    freeze=struct('status','ALL_BASELINES_AND_SUPPORTS_FROZEN_TEST_CLOSED', ...
        'PreregistrationSHA256',char(prereg_hash),'BaselineDefinitionsSHA256',char(definitions_hash), ...
        'PreTestAmendment01SHA256',char(amendment_hash),'PreTestAmendment02SHA256',char(amendment_02_hash), ...
        'PreTestAmendment03SHA256',char(amendment_03_hash), ...
        'InputHashesSHA256',r22a_io('hash',fullfile(output,'r22h_input_hashes.csv')), ...
        'ImplementationHashesSHA256',r22a_io('hash',fullfile(output,'r22h_implementation_hashes.csv')), ...
        'SupportHashesSHA256',r22a_io('hash',fullfile(output,'r22h_support_hashes.csv')), ...
        'SupportOverlapSHA256',r22a_io('hash',fullfile(output,'r22h_support_overlap.csv')), ...
        'RunnerSHA256',char(r22g_math('sha_file',fullfile(audit,'run_r22h.m'))), ...
        'TEST_ACCESS_COUNT',0,'TEST_METRICS_ACCESSED',false,'NewChannelsGenerated',false);
    r22a_io('json',fullfile(output,'r22h_test_freeze.json'),freeze);
    status=struct('status','PREPARED_TEST_CLOSED','TEST_ACCESS_COUNT',0,'TEST_METRICS_ACCESSED',false, ...
        'NewChannelsGenerated',false,'PreregistrationSHA256',char(prereg_hash), ...
        'DefinitionCheck',definition_check,'R22GIndependentSupportSHA256',char(support_hash(5)));
    r22a_io('json',fullfile(output,'r22h_execution_status.json'),status);
    receipt=struct('status','PREPARED_TEST_CLOSED','TRAINChannelsVerifiedAndUsed',200, ...
        'VALChannelsVerified',100,'TESTChannelsVerifiedWithoutMetricAccess',1000, ...
        'R22DArtifactsVerifiedExceptDocumentedPostprocessor',true,'R22GArtifactManifestVerified',true, ...
        'IndependentSupportMatchesR22G',true,'HaarSupportMatchesR22G',true, ...
        'PreTestAmendments',3,'DefinitionCheck',definition_check,'TEST_OPENED',false);
    r22a_io('json',fullfile(output,'r22h_prepare_receipt.json'),receipt);
    result=receipt;
end

function result=refresh_pretest_freeze(root,audit,output,r22d_results,r22g_results,p)
    status=jsondecode(fileread(fullfile(output,'r22h_execution_status.json')));
    old_freeze=jsondecode(fileread(fullfile(output,'r22h_test_freeze.json')));
    amendment=jsondecode(fileread(fullfile(audit,'r22h_preregistration_amendment_03.json')));
    assert(strcmp(status.status,'PREPARED_TEST_CLOSED')&&status.TEST_ACCESS_COUNT==0 ...
        &&~status.TEST_METRICS_ACCESSED&&strcmp(status.PreregistrationSHA256,amendment.previous_preregistration_sha256), ...
        'R22H:RefreshState','Only the untouched pre-TEST freeze can be refreshed.');
    assert(strcmp(old_freeze.SupportHashesSHA256,r22a_io('hash',fullfile(output,'r22h_support_hashes.csv'))) ...
        &&strcmp(old_freeze.SupportOverlapSHA256,r22a_io('hash',fullfile(output,'r22h_support_overlap.csv'))), ...
        'R22H:RefreshSupports','Frozen supports or overlap changed; refusing refresh.');
    assert(strcmp(amendment.amended_preregistration_sha256, ...
        r22g_math('sha_file',fullfile(audit,'r22h_preregistration.json'))));
    assert(strcmp(amendment.previous_preregistration_sha256,old_freeze.PreregistrationSHA256));
    verify_frozen_inputs(root,r22d_results,r22g_results,p);
    frozen=load(fullfile(output,'r22h_frozen_transforms.mat'),'U','masks');
    supports=r22g_math('read',fullfile(output,'r22h_support_hashes.csv'));
    for method=1:5
        assert(nnz(frozen.masks{method})==p.support_protocol.K ...
            &&string(r22g_math('support_hash',frozen.masks{method}))==supports.MaskSHA256(method));
    end
    assert(~isfile(fullfile(output,'r22h_reproduction_gate_partial.csv')) ...
        &&~isfile(fullfile(output,'r22h_test_per_channel.csv')), ...
        'R22H:RefreshAfterTest','TEST artifacts exist; refusing pre-TEST refresh.');
    copyfile(fullfile(audit,'r22h_preregistration.json'),fullfile(output,'r22h_preregistration.json'),'f');
    copyfile(fullfile(audit,'r22h_preregistration_amendment_03.json'), ...
        fullfile(output,'r22h_preregistration_amendment_03.json'),'f');
    implementation=make_implementation_hashes(root,audit);
    writetable(implementation,fullfile(output,'r22h_implementation_hashes.csv'));
    freeze=old_freeze;freeze.PreregistrationSHA256=amendment.amended_preregistration_sha256;
    freeze.PreTestAmendment03SHA256=r22g_math('sha_file', ...
        fullfile(audit,'r22h_preregistration_amendment_03.json'));
    freeze.ImplementationHashesSHA256=r22a_io('hash',fullfile(output,'r22h_implementation_hashes.csv'));
    freeze.RunnerSHA256=r22g_math('sha_file',fullfile(audit,'run_r22h.m'));
    r22a_io('json',fullfile(output,'r22h_test_freeze.json'),freeze);
    status.PreregistrationSHA256=amendment.amended_preregistration_sha256;status.PreTestAmendments=3;
    status.TEST_ACCESS_COUNT=0;status.TEST_METRICS_ACCESSED=false;
    r22a_io('json',fullfile(output,'r22h_execution_status.json'),status);
    receipt=jsondecode(fileread(fullfile(output,'r22h_prepare_receipt.json')));
    receipt.PreTestAmendments=3;receipt.FreezeRefreshedAfterValidationRepair=true;
    r22a_io('json',fullfile(output,'r22h_prepare_receipt.json'),receipt);
    result=struct('status','PREPARED_TEST_CLOSED','SupportsRecomputed',false, ...
        'PreTestAmendments',3,'TEST_ACCESS_COUNT',0,'TEST_METRICS_ACCESSED',false);
end

function result=evaluate_r22h(root,audit,output,r22d_results,r22g_results,p)
    status=jsondecode(fileread(fullfile(output,'r22h_execution_status.json')));
    freeze=jsondecode(fileread(fullfile(output,'r22h_test_freeze.json')));
    assert(strcmp(status.status,'PREPARED_TEST_CLOSED')&&status.TEST_ACCESS_COUNT==0 ...
        &&~status.TEST_METRICS_ACCESSED,'R22H:TestState','TEST is not closed in the expected frozen state.');
    assert(strcmp(freeze.PreregistrationSHA256,r22g_math('sha_file',fullfile(audit,'r22h_preregistration.json'))));
    assert(strcmp(freeze.BaselineDefinitionsSHA256,r22g_math('sha_file',fullfile(audit,'r22h_baseline_definitions.md'))));
    assert(strcmp(freeze.PreTestAmendment01SHA256,r22g_math('sha_file', ...
        fullfile(audit,'r22h_preregistration_amendment_01.json'))));
    assert(strcmp(freeze.PreTestAmendment02SHA256,r22g_math('sha_file', ...
        fullfile(audit,'r22h_preregistration_amendment_02.json'))));
    assert(strcmp(freeze.PreTestAmendment03SHA256,r22g_math('sha_file', ...
        fullfile(audit,'r22h_preregistration_amendment_03.json'))));
    assert(strcmp(freeze.RunnerSHA256,r22g_math('sha_file',fullfile(audit,'run_r22h.m'))));
    assert(strcmp(freeze.SupportHashesSHA256,r22a_io('hash',fullfile(output,'r22h_support_hashes.csv'))));
    assert(strcmp(freeze.InputHashesSHA256,r22a_io('hash',fullfile(output,'r22h_input_hashes.csv'))));
    assert(strcmp(freeze.ImplementationHashesSHA256,r22a_io('hash',fullfile(output,'r22h_implementation_hashes.csv'))));
    check_implementation_hashes(fullfile(output,'r22h_implementation_hashes.csv'));
    hashes=verify_frozen_inputs(root,r22d_results,r22g_results,p);
    [names,~,U,~,~,masks]=method_specs(r22g_results);
    stored=load(fullfile(output,'r22h_frozen_transforms.mat'),'U','masks');
    for method=1:5
        assert(norm(U{method}-stored.U{method},'fro')<1e-14, ...
            'R22H:FrozenTransformChanged','Frozen transform changed for %s.',names(method));
    end
    assert(isequal(masks{4},stored.masks{4})&&isequal(masks{5},stored.masks{5}), ...
        'R22H:FrozenReferenceSupportChanged','Frozen Haar/Independent reference masks changed.');
    masks=stored.masks;
    supports=r22g_math('read',fullfile(output,'r22h_support_hashes.csv'));
    test_manifest=readtable(fullfile(r22d_results,'r22d_test_manifest.csv'),'TextType','string');
    artifact=r22g_math('read',fullfile(r22d_results,'artifact_manifest.csv'));
    status.status='TEST_OPENED_REPRODUCIBILITY_GATE';status.TEST_ACCESS_COUNT=1;status.TEST_METRICS_ACCESSED=true;
    status.R22DTestManifestSHA256=char(hashes.SHA256(hashes.Paths==string(fullfile(r22d_results,'r22d_test_manifest.csv'))));
    r22a_io('json',fullfile(output,'r22h_execution_status.json'),status);
    snr=p.physical_model.snr_db(:).';gate_methods=[5,4];gate_table=table();
    for channel=1:height(test_manifest)
        G=load_test_channel(r22d_results,test_manifest(channel,:),artifact);
        rows=evaluate_channel(G,U,masks,gate_methods,names,test_manifest(channel,:),snr);
        if isempty(gate_table),gate_table=rows;else,gate_table=[gate_table;rows];end %#ok<AGROW>
        if mod(channel,50)==0||channel==height(test_manifest)
            writetable(gate_table,fullfile(output,'r22h_reproduction_gate_partial.csv'));
            status.GateChannelsCompleted=channel;r22a_io('json',fullfile(output,'r22h_execution_status.json'),status);
            fprintf('R22H reproduction gate %d/1000 channels\n',channel);
        end
    end
    gate_table=sortrows(gate_table,{'VelocityIndex','ChannelID','Representation'});
    writetable(gate_table,fullfile(output,'r22h_reproduction_gate.csv'));
    stored_summary=r22g_math('read',fullfile(r22g_results,'r22g_test_summary.csv'));
    expected_names=["I-R";"H-H"];expected_values=[p.reproducibility_gates.independent_exact_expected_excess_mse; ...
        p.reproducibility_gates.haar_exact_expected_excess_mse];
    gate_rows=cell(2,5);gate_pass=true;
    for index=1:2
        name=names(gate_methods(index));subset=gate_table(gate_table.Representation==name,:);
        fresh=mean(subset.ExcessAggregate);prior=stored_summary.MeanExcessMSE(stored_summary.Method==expected_names(index));
        assert(isscalar(prior)&&abs(prior-expected_values(index))<1e-14);
        difference=abs(fresh-expected_values(index));passed=difference<=p.reproducibility_gates.mean_absolute_tolerance;
        gate_pass=gate_pass&&passed;
        gate_rows(index,:)={name,fresh,prior,difference,passed};
    end
    gate_report=cell2table(gate_rows,'VariableNames',{'Representation','FreshTestMeanExcessMSE', ...
        'FrozenR22GMeanExcessMSE','AbsoluteDifference','Pass'});
    writetable(gate_report,fullfile(output,'r22h_reproducibility_gate.csv'));
    r22a_io('json',fullfile(output,'r22h_reproduction_gate_receipt.json'),struct( ...
        'Pass',gate_pass,'Tolerance',p.reproducibility_gates.mean_absolute_tolerance, ...
        'FrozenIndependentSupportSHA256',char(supports.MaskSHA256(supports.Representation==names(5))), ...
        'R22GIndependentSupportSHA256',p.reproducibility_gates.independent_support_sha256, ...
        'IndependentAndHaarOnlyEvaluatedBeforeGate',true));
    if ~gate_pass
        status.status='STOPPED_REPRODUCIBILITY_GATE';status.STOP_REASON='Fresh frozen-method TEST metrics did not reproduce R2.2G within tolerance.';
        status.ScientificClassification='RESULT_REQUIRES_AUDIT';r22a_io('json',fullfile(output,'r22h_execution_status.json'),status);
        write_stopped_report(output,gate_report);
        write_manifest(output,p,'RESULT_REQUIRES_AUDIT',hashes);
        result=struct('status','STOPPED','Classification','RESULT_REQUIRES_AUDIT','Gate',gate_report);return
    end
    method_rows=cell(1,5);method_rows{4}=gate_table(gate_table.Representation==names(4),:);
    method_rows{5}=gate_table(gate_table.Representation==names(5),:);
    other_methods=[1,2,3];other_table=table();
    for channel=1:height(test_manifest)
        G=load_test_channel(r22d_results,test_manifest(channel,:),artifact);
        rows=evaluate_channel(G,U,masks,other_methods,names,test_manifest(channel,:),snr);
        if isempty(other_table),other_table=rows;else,other_table=[other_table;rows];end %#ok<AGROW>
        if mod(channel,50)==0||channel==height(test_manifest)
            writetable(other_table,fullfile(output,'r22h_natural_baselines_partial.csv'));
            status.NaturalBaselineChannelsCompleted=channel;r22a_io('json',fullfile(output,'r22h_execution_status.json'),status);
            fprintf('R22H natural baselines %d/1000 channels\n',channel);
        end
    end
    for method=1:3,method_rows{method}=other_table(other_table.Representation==names(method),:);end
    test_table=vertcat(method_rows{:});test_table=sortrows(test_table,{'VelocityIndex','ChannelID','Representation'});
    assert(height(test_table)==5000&&all(test_table.K==p.support_protocol.K));
    writetable(test_table,fullfile(output,'r22h_test_per_channel.csv'));
    full_table=make_full_mse_table(test_table,snr);writetable(full_table,fullfile(output,'r22h_full_mse.csv'));
    [summary,per_channel,contrasts,classification]=summarize_results(test_table,names,p,supports);
    writetable(summary,fullfile(output,'r22h_test_summary.csv'));
    writetable(contrasts,fullfile(output,'r22h_paired_contrasts.csv'));
    diagnostics=make_diagnostics(test_table,supports);writetable(diagnostics,fullfile(output,'r22h_diagnostics.csv'));
    invariance=make_invariance(test_table,U,names);writetable(invariance,fullfile(output,'r22h_invariance.csv'));
    assert(max(invariance.MaxSingularValueRelativeDifference)<1e-10 ...
        &&max(invariance.MaxFullMSEAbsoluteDifference)<1e-10,'R22H:FullInvariance','Full-channel invariance failed.');
    assert(all(isfinite(test_table.ExcessAggregate))&&min(test_table.ExcessMSE_SNR03)>-1e-10, ...
        'R22H:Metric','Nonfinite metric or sparse MSE below full MSE.');
    report=write_report(output,summary,contrasts,diagnostics,classification,gate_report);
    status.status='COMPLETED';status.ScientificClassification=classification;
    status.FullChannelInvariance='PASS';status.TESTRows=height(test_table);
    r22a_io('json',fullfile(output,'r22h_execution_status.json'),status);
    write_manifest(output,p,classification,hashes);
    result=struct('status','PASS','Summary',summary,'Contrasts',contrasts, ...
        'Classification',classification,'FullChannelInvariance','PASS','Report',report);
    assert(size(per_channel,1)==1000);
end

function hashes=verify_frozen_inputs(root,r22d_results,r22g_results,p)
    audit=fullfile(root,'pufb_otfs_v2','audit','r2_2h_natural_baselines');
    paths={fullfile(root,p.ensemble_hashes.TRAIN_manifest.path), ...
        fullfile(root,p.ensemble_hashes.VAL_manifest.path), ...
        fullfile(root,p.ensemble_hashes.TEST_manifest.path), ...
        fullfile(root,p.ensemble_hashes.R22D_artifact_manifest.path), ...
        fullfile(root,p.ensemble_hashes.R22G_freeze.path), ...
        fullfile(root,p.ensemble_hashes.R22G_support_hashes.path), ...
        fullfile(root,p.ensemble_hashes.R22G_test_per_channel.path), ...
        fullfile(root,p.ensemble_hashes.R22G_artifact_manifest.path)};
    expected=[string(p.ensemble_hashes.TRAIN_manifest.sha256);string(p.ensemble_hashes.VAL_manifest.sha256); ...
        string(p.ensemble_hashes.TEST_manifest.sha256);string(p.ensemble_hashes.R22D_artifact_manifest.sha256); ...
        string(p.ensemble_hashes.R22G_freeze.sha256);string(p.ensemble_hashes.R22G_support_hashes.sha256); ...
        string(p.ensemble_hashes.R22G_test_per_channel.sha256);string(p.ensemble_hashes.R22G_artifact_manifest.sha256)];
    actual=strings(numel(paths),1);
    for index=1:numel(paths)
        actual(index)=string(r22g_math('sha_file',paths{index}));
        assert(actual(index)==expected(index),'R22H:InputHash','Frozen input hash mismatch: %s',paths{index});
    end
    amendment=jsondecode(fileread(fullfile(audit,'r22h_preregistration_amendment_01.json')));
    amendment_02=jsondecode(fileread(fullfile(audit,'r22h_preregistration_amendment_02.json')));
    amendment_03=jsondecode(fileread(fullfile(audit,'r22h_preregistration_amendment_03.json')));
    amendment_hash=string(r22g_math('sha_file',fullfile(audit,'r22h_preregistration.json')));
    assert(string(amendment.amended_preregistration_sha256)==string(amendment_02.previous_preregistration_sha256) ...
        &&string(amendment_02.amended_preregistration_sha256)==string(amendment_03.previous_preregistration_sha256) ...
        &&amendment_hash==string(amendment_03.amended_preregistration_sha256));
    assert(string(amendment_02.previous_baseline_definitions_sha256) ...
        =="824620ea7d0d5b5db842859770f83aa51063fecc4024a47e0a6405b2583d5dd5" ...
        &&string(amendment_02.amended_baseline_definitions_sha256) ...
        ==string(r22g_math('sha_file',fullfile(audit,'r22h_baseline_definitions.md'))));
    r22d_artifacts=r22g_math('read',fullfile(r22d_results,'artifact_manifest.csv'));
    exception_path=string(fullfile(root,amendment.exception.path));exception_seen=false;
    for index=1:height(r22d_artifacts)
        path=string(r22d_artifacts.AbsolutePath(index));expected_hash=string(r22d_artifacts.SHA256(index));
        actual_hash=string(r22g_math('sha_file',char(path)));
        if path==exception_path
            assert(expected_hash==string(amendment.exception.frozen_sha256) ...
                &&actual_hash==string(amendment.exception.observed_sha256), ...
                'R22H:AmendmentScope','R2.2D exception no longer matches its PRE-TEST amendment.');
            exception_seen=true;paths{end+1}=char(path);actual(end+1,1)=actual_hash; %#ok<AGROW>
        else
            assert(actual_hash==expected_hash,'R22H:ArtifactHash','R2.2D artifact hash mismatch: %s',path);
        end
    end
    assert(exception_seen,'R22H:AmendmentScope','Documented R2.2D postprocessor exception was not found.');
    r22a_io('verify',fullfile(r22g_results,'artifact_manifest.csv'));
    r22a_io('verify',fullfile(r22d_results,'r22d_test_input_hashes.csv'));
    sources=r22g_math('read',fullfile(r22g_results,'r22g_source_hashes.csv'));
    for index=1:height(sources)
        actual_hash=string(r22g_math('sha_file',char(sources.AbsolutePath(index))));
        assert(actual_hash==sources.FrozenSHA256(index)&&actual_hash==sources.ActualSHA256(index), ...
            'R22H:CriticalSourceHash','Critical source hash mismatch: %s',sources.AbsolutePath(index));
    end
    hashes=table(string(paths(:)),actual,'VariableNames',{'Paths','SHA256'});
end

function validate_manifest(manifest,count,seed_base,per_velocity)
    velocity=repelem((1:5)',per_velocity);channel=repmat((1:per_velocity)',5,1);
    assert(height(manifest)==count&&all(manifest.ChannelGenerated));
    assert(isequal(manifest.VelocityIndex,velocity)&&isequal(manifest.ChannelID,channel));
    assert(isequal(manifest.Seed,seed_base+1000*velocity+channel));
end

function check=check_baseline_definitions()
    M=24;N=12;F_M=dftmtx(M)/sqrt(M);F_N=dftmtx(N)/sqrt(N);
    X=reshape((1:M*N)+1i*cos(1:M*N),M,N);TF=F_M*X*F_N';S=F_M'*TF;
    U_fourier=kron(F_N',eye(M));equivalence=norm(S(:)-U_fourier*X(:))/norm(S(:));
    fourier_unitarity=norm(U_fourier'*U_fourier-eye(M*N),'fro');
    [~,~,haar]=fir_mallat_r22a(M,pi/4);[~,~,d4]=fir_mallat_r22a(M,pi/3);
    expected_haar=[1,1,0,0]/sqrt(2);sqrt_three=sqrt(3);
    expected_d4=[1+sqrt_three,3+sqrt_three,3-sqrt_three,1-sqrt_three]/(4*sqrt(2));
    expected_d4_high=(-1).^(0:3).*fliplr(expected_d4);
    check=struct('FourierExplicitRelativeError',equivalence,'FourierUnitarityResidual',fourier_unitarity, ...
        'HaarLowpassError',norm(haar.h-expected_haar),'HaarHighpassError', ...
        norm(haar.g-[0,0,1,-1]/sqrt(2)),'D4LowpassError',norm(d4.h-expected_d4), ...
        'D4HighpassError',norm(d4.g-expected_d4_high),'HaarOrthogonalityResidual',haar.OrthogonalityResidual, ...
        'D4OrthogonalityResidual',d4.OrthogonalityResidual);
    assert(equivalence<1e-12&&fourier_unitarity<1e-12&&check.HaarLowpassError<1e-14 ...
        &&check.HaarHighpassError<1e-14&&check.D4LowpassError<1e-14&&check.D4HighpassError<1e-14 ...
        &&haar.OrthogonalityResidual<1e-12&&d4.OrthogonalityResidual<1e-12, ...
        'R22H:BaselineDefinition','Fourier, Haar or D4 definition gate failed.');
end

function [names,definitions,U,alpha_nu,alpha_tau,frozen_masks]=method_specs(r22g_results)
    M=24;N=12;F_N=dftmtx(N)/sqrt(N);
    [Pnu_haar,~,~]=fir_mallat_r22a(N,pi/4);[Ptau_haar,~,~]=fir_mallat_r22a(M,pi/4);
    [Pnu_d4,~,~]=fir_mallat_r22a(N,pi/3);[Ptau_d4,~,~]=fir_mallat_r22a(M,pi/3);
    frozen=load(fullfile(r22g_results,'r22g_selected_methods.mat'),'method_cells');
    method_cells=frozen.method_cells;independent=method_cells{4};haar_frozen=method_cells{5};
    assert(string(independent.Name)=="I-R"&&string(haar_frozen.Name)=="H-H");
    U={eye(M*N),kron(F_N',eye(M)),kron(Pnu_d4,Ptau_d4), ...
        kron(Pnu_haar,Ptau_haar),independent.U};
    definitions=["U=I_288";"U=kron(F_N^H,I_M), audited conventional OTFS DD-to-time map"; ...
        "kron(P_nu(pi/3),P_tau(pi/3)), periodic D4 partial Mallat"; ...
        "kron(P_nu(pi/4),P_tau(pi/4)), periodic Haar partial Mallat"; ...
        "frozen R2.2G I-R, (alpha_nu,alpha_tau)=(0,pi/4)"];
    names=["Identity";"Fourier OTFS";"D4xD4";"HaarxHaar";"Independent PUFB"];
    alpha_nu=[NaN;NaN;pi/3;pi/4;0];alpha_tau=[NaN;NaN;pi/3;pi/4;pi/4];
    frozen_masks={[],[],[],haar_frozen.Mask,independent.Mask};
    current=r22d_model('transform',0,pi/4);
    assert(norm(current.U-independent.U,'fro')<1e-13);
    assert(norm(kron(Pnu_haar,Ptau_haar)-haar_frozen.U,'fro')<1e-13);
    for index=1:5,assert(norm(U{index}'*U{index}-eye(M*N),'fro')<1e-11);end
end

function bank=load_partition(r22d_results,manifest,split)
    bank=complex(zeros(288,288,height(manifest)));
    for index=1:height(manifest)
        if strcmp(split,'train')
            path=fullfile(r22d_results,'train_channels',sprintf('v%d_c%03d.mat', ...
                manifest.VelocityIndex(index),manifest.ChannelID(index)));
        else
            error('R22H:Split','Unsupported TRAIN loading split.');
        end
        saved=load(path,'G');bank(:,:,index)=saved.G;
    end
end

function G=load_test_channel(r22d_results,row,artifact)
    channel=(row.VelocityIndex-1)*200+row.ChannelID;
    path=fullfile(r22d_results,'test_channels',sprintf('channel_%04d.mat',channel));
    record=artifact(artifact.AbsolutePath==string(path),:);
    assert(height(record)==1&&strcmp(r22a_io('hash',path),record.SHA256), ...
        'R22H:TestChannelHash','TEST channel hash mismatch: %s',path);
    saved=load(path,'G');G=saved.G;
end

function rows=evaluate_channel(G,U,masks,methods,names,manifest_row,snr)
    n=size(G,1);reference_singular=svd(G);reference_full=zeros(1,numel(snr));
    for point=1:numel(snr)
        variance=10^(-snr(point)/10);reference_full(point)=mean(variance./(reference_singular.^2+variance));
    end
    rows=[];
    for method=methods
        H=U{method}'*G*U{method};singular=svd(H);singular_difference=norm(singular-reference_singular)/max(norm(reference_singular),eps);
        full=zeros(1,numel(snr));sparse=zeros(1,numel(snr));excess=zeros(1,numel(snr));
        Hs=H.*masks{method};gram=Hs'*Hs;identity=eye(n);
        for point=1:numel(snr)
            variance=10^(-snr(point)/10);factor=chol(gram+variance*identity);
            W=factor\(factor'\Hs');residual=W*H-identity;
            sparse(point)=(norm(residual,'fro')^2+variance*norm(W,'fro')^2)/n;
            full(point)=mean(variance./(singular.^2+variance));excess(point)=sparse(point)-full(point);
        end
        full_difference=max(abs(full-reference_full));energy=abs(H).^2;total=sum(energy,'all');probability=energy(:)/total;
        [ordered,~]=sort(probability,'descend');K99=find(cumsum(ordered)>=.99,1);
        omitted=sum(energy(~masks{method}),'all')/total;
        sparse_aggregate=regime_mean(sparse,snr);excess_aggregate=regime_mean(excess,snr);
        assert(all(isfinite([sparse,full,excess]))&&min(excess)>-1e-10, ...
            'R22H:Metric','Nonfinite or invalid MSE for %s.',names(method));
        metadata={names(method),manifest_row.VelocityIndex,manifest_row.ChannelID,manifest_row.Seed, ...
            sparse_aggregate,excess_aggregate,K99,omitted,1-omitted,33178, ...
            string(r22g_math('support_hash',masks{method})),singular_difference,full_difference};
        row=cell2table(metadata,'VariableNames',{'Representation','VelocityIndex','ChannelID','Seed', ...
            'SparseAggregate','ExcessAggregate','K99','OmittedEnergyFraction','RetainedEnergyFraction','K', ...
            'SupportSHA256','SingularValueRelativeDifference','MaxFullMSEAbsoluteDifference'});
        for point=1:numel(snr)
            row.(sprintf('SparseMSE_SNR%02d',snr(point)))=sparse(point);
            row.(sprintf('FullMSE_SNR%02d',snr(point)))=full(point);
            row.(sprintf('ExcessMSE_SNR%02d',snr(point)))=excess(point);
        end
        if isempty(rows),rows=row;else,rows=[rows;row];end %#ok<AGROW>
    end
end

function value=regime_mean(metrics,snr)
    low=ismember(snr,[3,6,9]);mid=ismember(snr,[12,15,18]);high=ismember(snr,[21,24]);
    value=mean([mean(metrics(low)),mean(metrics(mid)),mean(metrics(high))]);
end

function full_table=make_full_mse_table(test_table,snr)
    rows=cell(height(test_table)*numel(snr),6);cursor=0;
    for index=1:height(test_table)
        for point=1:numel(snr)
            cursor=cursor+1;column=sprintf('FullMSE_SNR%02d',snr(point));
            rows(cursor,:)={test_table.Representation(index),test_table.VelocityIndex(index), ...
                test_table.ChannelID(index),test_table.Seed(index),snr(point),test_table.(column)(index)};
        end
    end
    full_table=cell2table(rows,'VariableNames',{'Representation','VelocityIndex','ChannelID','Seed','SNRdB','FullMSE'});
end

function [summary,per_channel,contrasts,classification]=summarize_results(test_table,names,p,supports)
    per_channel=zeros(1000,5);sparse_channel=per_channel;full_channel=per_channel;
    for method=1:5
        subset=sortrows(test_table(test_table.Representation==names(method),:),{'VelocityIndex','ChannelID'});
        assert(height(subset)==1000);per_channel(:,method)=subset.ExcessAggregate;
        sparse_channel(:,method)=subset.SparseAggregate;
        full_columns=arrayfun(@(v)sprintf('FullMSE_SNR%02d',v),p.physical_model.snr_db,'UniformOutput',false);
        full_channel(:,method)=arrayfun(@(row)regime_mean(subset{row,full_columns},p.physical_model.snr_db),1:1000).';
    end
    [means,lower,upper]=r22g_math('bootstrap',per_channel,p.bootstrap.seed,p.bootstrap.replicates);
    records=cell(5,8);
    for method=1:5
        item=supports(supports.Representation==names(method),:);
        records(method,:)={names(method),means(method),lower(method),upper(method),mean(sparse_channel(:,method)), ...
            mean(full_channel(:,method)),p.support_protocol.K,item.MaskSHA256};
    end
    summary=cell2table(records,'VariableNames',{'Representation','MeanExcessMSE','CI95Lower','CI95Upper', ...
        'MeanSparseMSE','MeanFullMSE','K','SupportSHA256'});
    contrast_pairs=[5,1;5,2;5,3;5,4];labels=["INDEPENDENT_MINUS_IDENTITY";"INDEPENDENT_MINUS_FOURIER_OTFS"; ...
        "INDEPENDENT_MINUS_D4";"INDEPENDENT_MINUS_HAAR"];
    differences=zeros(1000,4);
    for contrast=1:4,differences(:,contrast)=per_channel(:,contrast_pairs(contrast,1))-per_channel(:,contrast_pairs(contrast,2));end
    [difference_mean,difference_lower,difference_upper]=r22g_math('bootstrap',differences,p.bootstrap.seed,p.bootstrap.replicates);
    rows=cell(4,6);
    for contrast=1:4
        baseline=contrast_pairs(contrast,2);relative=-difference_mean(contrast)/means(baseline);
        rows(contrast,:)={labels(contrast),difference_mean(contrast),difference_lower(contrast), ...
            difference_upper(contrast),relative,names(baseline)};
    end
    contrasts=cell2table(rows,'VariableNames',{'Contrast','MeanPairedDifference','CI95Lower','CI95Upper','RelativeReduction','Baseline'});
    classification=classify_result(difference_lower,difference_upper);
end

function classification=classify_result(lower,upper)
    if any(lower>0)
        classification="NATURAL_BASELINE_OUTPERFORMS_INDEPENDENT";
    elseif all(upper<0)
        classification="INDEPENDENT_BEATS_ALL_NATURAL_BASELINES";
    elseif upper(3)<0&&upper(4)<0&&~(upper(1)<0&&upper(2)<0)
        classification="INDEPENDENT_BEATS_FIXED_WAVELET_BASELINES_ONLY";
    else
        classification="INDEPENDENT_COMPARABLE_TO_NATURAL_BASELINES";
    end
end

function diagnostics=make_diagnostics(test_table,supports)
    names=supports.Representation;rows=cell(5,7);
    for method=1:5
        subset=test_table(test_table.Representation==names(method),:);
        rows(method,:)={names(method),mean(subset.K99),mean(subset.OmittedEnergyFraction), ...
            mean(subset.RetainedEnergyFraction),supports.TrainRetainedEnergyFraction(method), ...
            supports.K(method),supports.MaskSHA256(method)};
    end
    diagnostics=cell2table(rows,'VariableNames',{'Representation','MeanK99','MeanOmittedEnergyFraction', ...
        'MeanRetainedEnergyFraction','TrainRetainedEnergyFraction','K','SupportSHA256'});
end

function invariance=make_invariance(test_table,U,names)
    rows=cell(numel(names),5);
    for method=1:numel(names)
        subset=test_table(test_table.Representation==names(method),:);
        rows(method,:)={names(method),height(subset),max(subset.SingularValueRelativeDifference), ...
            max(subset.MaxFullMSEAbsoluteDifference),norm(U{method}'*U{method}-eye(288),'fro')};
    end
    invariance=cell2table(rows,'VariableNames',{'Representation','ChannelCount', ...
        'MaxSingularValueRelativeDifference','MaxFullMSEAbsoluteDifference','UnitarityResidual'});
end

function report=write_report(output,summary,contrasts,diagnostics,classification,gate)
    lines=["# R2.2H Natural Baselines";""; ...
        "R22H_STATUS = PASS. All five sparse representations were evaluated on the same frozen R2.2D TEST channels with K=33178 and lambda=sigma^2. No new channels, BER, eta, K sweep, or manuscript changes.";""; ...
        "## Sparse Receiver Results";""; ...
        "| Representation | Definition | Excess MSE | 95% CI | Mean sparse MSE | Relative excess vs Independent |"; ...
        "|---|---|---:|---:|---:|---:|"];
    independent=summary.MeanExcessMSE(summary.Representation=="Independent PUFB");
    for method=1:height(summary)
        relative=(summary.MeanExcessMSE(method)-independent)/independent;
        lines(end+1)=sprintf('| %s | %s | %.12g | [%.12g, %.12g] | %.12g | %.6g |', ...
            char(summary.Representation(method)),char(definition_for(summary.Representation(method))), ...
            summary.MeanExcessMSE(method),summary.CI95Lower(method),summary.CI95Upper(method), ...
            summary.MeanSparseMSE(method),relative); %#ok<AGROW>
    end
    lines=[lines;"";"## Paired Contrasts";""; ...
        "Sign convention: Independent minus baseline; negative favors Independent. Relative reduction is `-contrast / baseline mean excess MSE`.";""; ...
        "| Contrast | Mean paired difference | 95% paired bootstrap CI | Relative reduction |"; ...
        "|---|---:|---:|---:|"];
    for index=1:height(contrasts)
        lines(end+1)=sprintf('| %s | %.12g | [%.12g, %.12g] | %.6g |',char(contrasts.Contrast(index)), ...
            contrasts.MeanPairedDifference(index),contrasts.CI95Lower(index),contrasts.CI95Upper(index), ...
            contrasts.RelativeReduction(index)); %#ok<AGROW>
    end
    lines=[lines;"";"## Full-Channel Reference and Diagnostics";""; ...
        "Full-channel LMMSE is reported as absolute MSE in `r22h_full_mse.csv`; it is not a sparse method and has no excess-MSE contrast against itself."; ...
        "Full-channel singular-value and MSE invariance: PASS at 1e-10 for every representation and TEST channel.";""; ...
        "| Representation | K99 | Omitted energy | Retained energy | Support SHA256 |";"|---|---:|---:|---:|---|"];
    for index=1:height(diagnostics)
        lines(end+1)=sprintf('| %s | %.6g | %.9g | %.9g | `%s` |',char(diagnostics.Representation(index)), ...
            diagnostics.MeanK99(index),diagnostics.MeanOmittedEnergyFraction(index), ...
            diagnostics.MeanRetainedEnergyFraction(index),char(diagnostics.SupportSHA256(index))); %#ok<AGROW>
    end
    lines=[lines;"";sprintf('R2.2G Independent/Haar reproduction gate: PASS (tolerance 1e-10).'); ...
        sprintf('Scientific classification: **%s**.',classification);""; ...
        "No computational-speedup or universal-superiority claim is made. No next experiment was launched."];
    report=fullfile(output,'R22H_NATURAL_BASELINES_REPORT.md');r22a_io('text',report,join(lines,newline));
    writetable(gate,fullfile(output,'r22h_reproduction_gate.csv'));
end

function definition=definition_for(name)
    switch name
        case "Identity",definition="fixed I_288";
        case "Fourier OTFS",definition="fixed kron(F_Nu^H,I_M)";
        case "D4xD4",definition="fixed alpha=pi/3 on both axes";
        case "HaarxHaar",definition="fixed alpha=pi/4 on both axes";
        otherwise,definition="frozen (0,pi/4)";
    end
end

function write_stopped_report(output,gate)
    lines=["# R2.2H Natural Baselines";""; ...
        "R22H_STATUS = STOPPED. The independent/Haar reproducibility gate failed. Identity, Fourier OTFS, and D4 were not evaluated or compared.";""; ...
        sprintf('Classification: **RESULT_REQUIRES_AUDIT**.');""; ...
        "No interpretation or next experiment was produced."];
    r22a_io('text',fullfile(output,'R22H_NATURAL_BASELINES_REPORT.md'),join(lines,newline));
    writetable(gate,fullfile(output,'r22h_reproduction_gate.csv'));
end

function write_manifest(output,p,classification,hashes)
    files={'r22h_preregistration.json','r22h_preregistration_amendment_01.json', ...
        'r22h_preregistration_amendment_02.json', ...
        'r22h_baseline_definitions.md','r22h_input_hashes.csv','r22h_implementation_hashes.csv', ...
        'r22h_test_freeze.json','r22h_prepare_receipt.json','r22h_reproduction_gate_receipt.json', ...
        'r22h_support_hashes.csv','r22h_support_overlap.csv','r22h_diagnostics.csv','r22h_full_mse.csv', ...
        'r22h_invariance.csv','r22h_test_per_channel.csv','r22h_test_summary.csv', ...
        'r22h_paired_contrasts.csv','R22H_NATURAL_BASELINES_REPORT.md','r22h_reproduction_gate.csv'};
    artifacts=struct('File',{},'SHA256',{});
    for index=1:numel(files)
        path=fullfile(output,files{index});
        if isfile(path),artifacts(end+1)=struct('File',files{index}, ...
                'SHA256',char(r22g_math('sha_file',path)));end %#ok<AGROW>
    end
    final_status="COMPLETED";if classification=="RESULT_REQUIRES_AUDIT",final_status="STOPPED";end
    manifest=struct('status',char(final_status),'ScientificClassification',char(classification), ...
        'PreregistrationSHA256',char(r22g_math('sha_file',fullfile(output,'r22h_preregistration.json'))), ...
        'InputHashes',table_to_struct(hashes),'Artifacts',artifacts,'TESTChannels',1000,'NewChannelsGenerated',false, ...
        'TESTAccessCount',1,'ManuscriptModified',false,'K',p.support_protocol.K);
    r22a_io('json',fullfile(output,'r22h_manifest.json'),manifest);
end

function output=table_to_struct(input)
    paths=cellstr(input.Paths);hashes=cellstr(input.SHA256);
    output=struct('Path',paths,'SHA256',hashes);
end

function implementation=make_implementation_hashes(root,audit)
    pkg=fileparts(fileparts(mfilename('fullpath')));
    paths=string({fullfile(pkg,'07_natural_baselines','run_r22h.m'); ...
        fullfile(pkg,'02_shared_lib','fir_mallat_r22a.m'); ...
        fullfile(pkg,'03_axis_specific_pufb','r22d_model.m'); ...
        fullfile(pkg,'01_channel_generator','generate_channel_r20a.m'); ...
        fullfile(root,'pufb_otfs_v2','audit','r1_0r_wotfs_reproduction','waveform_operators_r10r.m'); ...
        fullfile(root,'pufb_otfs_v2','audit','r1_0r_wotfs_reproduction','modulate_explicit_2d_r10r.m'); ...
        fullfile(root,'pufb_otfs_v2','audit','r1_0r_wotfs_reproduction','validate_kron_equivalence_r10r.m'); ...
        fullfile(root,'pufb_otfs_v2','audit','r2_2g_factorial_ablation','r22g_math.m')});
    hashes=strings(size(paths));
    for index=1:numel(paths),hashes(index)=string(r22g_math('sha_file',char(paths(index))));end
    implementation=table(paths,hashes,'VariableNames',{'AbsolutePath','SHA256'});
end

function check_implementation_hashes(path)
    implementation=r22g_math('read',path);
    for index=1:height(implementation)
        actual=string(r22g_math('sha_file',char(implementation.AbsolutePath(index))));
        assert(actual==implementation.SHA256(index),'R22H:ImplementationHash', ...
            'Implementation changed after TEST freeze: %s',implementation.AbsolutePath(index));
    end
end