function result=run_r22dhma(root,stage)
    pkg=fileparts(mfilename('fullpath'));pkg=fileparts(pkg);
    audit=fullfile(pkg,'05_highsnr_mismatch_audit');
    output=fullfile(root,'pufb_otfs_v2','results','r2_2d_highsnr_mismatch_audit');
    source=fullfile(root,'pufb_otfs_v2','results','r2_2d_ber_ext');
    addpath(audit,fullfile(pkg,'04_ber_confirmation'),fullfile(pkg,'02_shared_lib'));
    if ~isfolder(output),mkdir(output);end
    if strcmp(stage,'PREPARE'),result=prepare(audit,output,source);return;end
    assert(strcmp(stage,'RUN'),'R22DHMA:Stage','Use PREPARE or RUN.');
    r22a_io('verify',fullfile(output,'r22dhma_input_hashes.csv'));r22a_io('verify',fullfile(output,'r22dhma_code_hashes.csv'));
    p=jsondecode(fileread(fullfile(output,'r22dhma_preregistration.json')));p.snr_db=p.snr_db(:).';
    loaded=load(fullfile(source,'r22dbe_frozen_methods.mat'),'methods');methods=loaded.methods;
    manifest=r22dcb_math('read',fullfile(source,'r22dbe_test_manifest.csv'));
    records=cell(500,4);
    for channel=1:500
        persisted=load(fullfile(source,'phase_b_channels',sprintf('channel_%04d.mat',channel)),'G','seed');G=persisted.G;
        assert(persisted.seed==manifest.Seed(channel));local_channel=mod(channel-1,100)+1;
        for method=1:4
            path=fullfile(output,sprintf('metric_c%04d_m%d.mat',channel,method));
            if isfile(path),loaded=load(path,'record');record=loaded.record;
            else
                U=methods(method).U;H=U'*G*U;mask=methods(method).Mask;assert(nnz(mask)==p.K);
                record=r22dhma_kernel(H,mask,p.snr_db,p,local_channel<=10);save(path,'record','-v7.3');
            end
            records{channel,method}=record;
        end
        if mod(channel,10)==0||channel==500,fprintf('Mismatch exact diagnostic %d/500 existing channels\n',channel);end
    end
    tables=write_tables(records,methods,p,output);
    result=r22dhma_finish(methods,p,tables,source,output);
    r22a_io('verify',fullfile(output,'r22dhma_input_hashes.csv'));r22a_io('verify',fullfile(output,'r22dhma_code_hashes.csv'));
    r22a_io('json',fullfile(output,'r22dhma_execution_status.json'),struct('status','COMPLETED','NEW_CHANNELS_GENERATED',false, ...
        'BER_RERUN',false,'NEW_BITS_GENERATED',false,'NEW_NOISE_GENERATED',false));
    artifacts=r22a_io('manifest',output,audit);r22a_io('verify',fullfile(output,'artifact_manifest.csv'));
    fprintf('MISMATCH AUDIT COMPLETE %s; %d artifacts\n',result.PRIMARY_FAILURE_MECHANISM,height(artifacts));
end

function p=prepare(audit,output,source)
    assert(~isfile(fullfile(output,'r22dhma_preregistration.json')),'R22DHMA:Prepared','No overwriting preregistration.');
    for name={'artifact_manifest.csv','r22dbe_code_hashes.csv','r22dbe_input_hashes.csv','r22dbe_history_hashes.csv'},r22a_io('verify',fullfile(source,name{1}));end
    manifest=r22dcb_math('read',fullfile(source,'r22dbe_test_manifest.csv'));counts=r22dcb_math('read',fullfile(source,'r22dbe_ber_counts.csv'));
    channels=r22dcb_math('read',fullfile(source,'r22dbe_ber_per_channel.csv'));files=dir(fullfile(source,'phase_b_channels','channel_*.mat'));
    assert(numel(files)==500&&height(manifest)==500&&height(counts)==340&&height(channels)==34000);
    assert(isequal(manifest.VelocityIndex,repelem((1:5)',100))&&isequal(manifest.ChannelID,repmat((1:100)',5,1)));
    assert(isequal(unique(counts.SNRdB),(20:36)')&&all(counts.RetainedEntries==33178)&&all(channels.RetainedEntries==33178));
    loaded=load(fullfile(source,'r22dbe_frozen_methods.mat'),'methods');methods=loaded.methods;
    method_hashes=jsondecode(fileread(fullfile(source,'r22dbe_method_freeze.json')));
    support_hashes=jsondecode(fileread(fullfile(source,'r22dbe_support_freeze.json')));
    for method=1:4
        assert(nnz(methods(method).Mask)==33178&&strcmp(r22dcb_math('arrayhash',methods(method).U),method_hashes(method).U_SHA256) ...
            &&strcmp(r22dcb_math('arrayhash',uint8(methods(method).Mask(:))),support_hashes(method).Mask_SHA256));
    end
    for velocity=1:5
        assert(sum(manifest.VelocityIndex==velocity)==100);
        for snr=20:36
            T=counts(counts.VelocityIndex==velocity&counts.SNRdB==snr,:);assert(height(T)==4);
            for method=1:4
                Q=channels(channels.Method==T.Method(method)&channels.VelocityIndex==velocity&channels.SNRdB==snr,:);
                assert(height(Q)==100&&sum(Q.Errors)==T.Errors(method)&&sum(Q.Bits)==T.Bits(method));
            end
        end
    end
    frozen=r22dcb_math('read',fullfile(source,'artifact_manifest.csv'));
    writetable(frozen(:,{'AbsolutePath','SHA256'}),fullfile(output,'r22dhma_input_hashes.csv'));
    r22a_io('json',fullfile(output,'r22dhma_input_audit.json'),struct('HISTORICAL_INPUTS_VERIFIED',true,'Channels',500, ...
        'ChannelsPerVelocity',100,'BERCells',340,'BERChannelRows',34000,'SNRGrid',20:36,'K',33178, ...
        'CodeHashesPASS',true,'MethodHashesPASS',true,'SupportHashesPASS',true,'CountsPairingPASS',true));
    p=struct('campaign','R2.2D-HIGHSNR-MISMATCH-AUDIT','K',33178,'snr_db',20:36, ...
        'rank_relative_tolerance',1e-12,'identity_tolerance',1e-10,'monotonic_tolerance',1e-10, ...
        'backward_tolerance',1e-12,'svd_receiver_relative_tolerance',1e-8,'svd_mse_absolute_tolerance',1e-10, ...
        'condition_identity_relative_tolerance',1e-8,'reference_subset','First10channelIDs in eachvelocity; all4methods/all17SNRs', ...
        'mid_reference',[28,29,30],'high_region',[34,35,36],'bootstrap_seed',151100000,'bootstrap_replicates',10000, ...
        'bootstrap_rule','Paired physical channels; independently resample100 within velocity; equal5velocity means; same seed/resamples all contrasts', ...
        'mechanism_rule','Both primary methods require positive lowerCI mismatch and at leastone T_OMIT/WSE_F/Amplification_F contrast; full/solver/noise gates mustpass', ...
        'systematic_rule','At least95percentchannels perprimary have nondecreasing condA and W_F acrossall17SNRs within1e-10relative tolerance; high-mid increases supported; strongcondition ratio meanHIGH/meanMID>=2', ...
        'regularization_interaction_rule','Mismatchsupported; systematiccondA/Wgain; positive lowerCI condA andW_F; strongcondratio>=2; noiseHIGH-MID upperCI<0 forbothprimary', ...
        'axis_mechanism_rule','Smalleromitted:24dB eF upperCI<0; betterconditioning:36dB condHs andcondA upperCI<0; loweramplification:HIGH gain_F andWSE_F upperCI<0. Atleast2components MULTIFACTOR.', ...
        'dominance_rule','Largest absolute REG/OMIT/CROSS/NOISE aggregatedterm; fractions signed/notprobabilities; no singleterm if largestabsolute<.5sumabs acrossbothprimaryHIGH', ...
        'decision_precedence','Full/noise failure E; solverfailure D; supportedstronginteraction B; supportedmismatchandtrueHIGHincrease A; true sparsemonotonic C; no dominantterm F; otherwise UNRESOLVED', ...
        'correlation_rule','Pearson/Spearman17SNR-level descriptiveonly. ChannelSpearman36minus30changes; fivevelocity-specific and pooled equal-size strata; no causal pvalues.', ...
        'NEW_CHANNELS_GENERATED',false,'BER_RERUN',false,'NEW_TRANSFORM_TRAINED',false,'ALPHA_CHANGED',false, ...
        'SUPPORT_REBUILT',false,'K_CHANGED',false,'RECEIVER_CHANGED',false,'FILTER_LENGTH_ESCALATED',false,'M_FB_CHANGED',false,'MANUSCRIPT_MODIFIED',false);
    r22a_io('json',fullfile(output,'r22dhma_preregistration.json'),p);
    code_hashes=r22dcb_math('read',fullfile(source,'r22dbe_code_hashes.csv'));
    location=find(endsWith(code_hashes.AbsolutePath,'r22dcb_ber.m'));assert(isscalar(location));
    ber_source=fileread(char(code_hashes.AbsolutePath(location)));sigma2=10.^(-(20:36)/10);
    mapping=contains(ber_source,'variance=10^(-snr/10)')&&contains(ber_source,"r22dcb_math('receiver',Hs,variance)") ...
        &&contains(ber_source,'physical_noise=sqrt(variance/2)')&&all(diff(sigma2)<0);
    assert(mapping,'R22DHMA:Noise','Frozen BER source convention failed.');
    r22a_io('json',fullfile(output,'r22dhma_snr_noise_audit.json'),struct('SNR_NOISE_MAPPING_PASS',mapping, ...
        'SNRGrid',20:36,'Sigma2',sigma2,'GenerationAndRegularizationSameVariance',true,'FrozenSourceSHA256',code_hashes.SHA256(location), ...
        'source_expressions',{{'variance=10^(-snr/10)','receiver(Hs,variance)','sqrt(variance/2)*(randn+1i*randn)'}},'CONVENTION_CHANGED',false));
    H=reshape(sin(1:144),12,12)+1i*reshape(cos(3:146),12,12)+eye(12);mask=mod(reshape(1:144,12,12),3)~=0;
    test=r22dhma_kernel(H,mask,20:36,p,true);assert(max(test.Values(:,19:20),[],'all')<p.identity_tolerance);
    r22dhma_finish('SELFTEST',p,[],source,output);
    files=dir(fullfile(audit,'*.m'));paths=cell(numel(files),1);
    for index=1:numel(files)
        paths{index}=fullfile(audit,files(index).name);d=checkcode(paths{index},'-id');
        for item=1:numel(d),fprintf('%s:%d %s\n',files(index).name,d(item).line,d(item).message);end
        assert(isempty(d),'R22DHMA:Checkcode','Source diagnostic before analysis.');
    end
    paths=[paths;{which('r22dcb_math');which('r22a_io');fullfile(output,'r22dhma_preregistration.json')}];
    hashes=strings(numel(paths),1);for index=1:numel(paths),hashes(index)=string(r22a_io('hash',paths{index}));end
    writetable(table(string(paths(:)),hashes,'VariableNames',{'AbsolutePath','SHA256'}),fullfile(output,'r22dhma_code_hashes.csv'));
    r22a_io('json',fullfile(output,'r22dhma_execution_status.json'),struct('status','INPUTS_VERIFIED_AND_RULES_FROZEN', ...
        'NEW_CHANNELS_GENERATED',false,'BER_RERUN',false));
end

function tables=write_tables(records,methods,p,output)
    names={'MSE_TRUE_SPARSE','MSE_MODEL_CONSISTENT','MISMATCH_PENALTY','T_REG','T_OMIT','T_CROSS','T_NOISE','MSE_FULL', ...
        'CondA','LambdaMin','LambdaMax','W_2','W_F','WSE_2','WSE_F','Amplification_F','Amplification_2','BackwardResidual', ...
        'ReconstructionError','PenaltyEqualityError','SVDReceiverRelativeError','SVDMSEDifference','ConditionIdentityError'};
    rows=cell(34000,27);static_rows=cell(2000,11);limit_rows=cell(2000,5);singular_rows=cell(2000*288,5);
    cursor=0;static_cursor=0;singular_cursor=0;
    for channel=1:500
        velocity=ceil(channel/100);local=mod(channel-1,100)+1;
        for method=1:4
            record=records{channel,method};name=methods(method).Name;static_cursor=static_cursor+1;
            static_rows(static_cursor,:)=[{name,velocity,local},num2cell(record.Static)];
            limit_rows(static_cursor,:)=[{name,velocity,local},num2cell(record.Limit)];
            for index=1:288
                singular_cursor=singular_cursor+1;singular_rows(singular_cursor,:)={name,velocity,local,index,record.Singular(index)};
            end
            for point=1:17,cursor=cursor+1;rows(cursor,:)=[{name,velocity,local,p.snr_db(point)},num2cell(record.Values(point,:))];end
        end
    end
    detailed=cell2table(rows,'VariableNames',[{'Method','VelocityIndex','ChannelID','SNRdB'},names]);
    static=cell2table(static_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SigmaMax','SigmaMin','CondHs','NumericalRank','RankTolerance','RelativeOmittedF','OmittedF','Omitted2'});
    singular=cell2table(singular_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SingularIndex','SingularValue'});
    limits=cell2table(limit_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','MSE_INF_APPROX','OMITTED_GAIN_INF'});
    writetable(singular,fullfile(output,'r22dhma_channel_singular_values.csv'));
    writetable(detailed(:,[1:12,23:24]),fullfile(output,'r22dhma_exact_mse_decomposition.csv'));
    writetable(detailed(:,[1:4,13:15]),fullfile(output,'r22dhma_conditioning.csv'));
    writetable(detailed(:,[1:4,16:21]),fullfile(output,'r22dhma_receiver_gain.csv'));
    writetable(detailed(:,[1:4,22]),fullfile(output,'r22dhma_solver_residuals.csv'));
    writetable(detailed(isfinite(detailed.SVDReceiverRelativeError),[1:4,25:27]),fullfile(output,'r22dhma_svd_reference_check.csv'));
    writetable(limits,fullfile(output,'r22dhma_highsnr_limit.csv'));writetable(static,fullfile(output,'r22dhma_channel_conditioning.csv'));
    save(fullfile(output,'r22dhma_diagnostic_bundle.mat'),'detailed','static','limits','-v7.3');
    tables=struct('Detailed',detailed,'Static',static,'Limits',limits);
end