function result=run_r22e(root,stage)
    pkg=fileparts(mfilename('fullpath'));pkg=fileparts(pkg);
    audit=fullfile(pkg,'06_receiver_floor');output=fullfile(root,'pufb_otfs_v2','results','r2_2e_receiver_floor');
    source=fullfile(root,'pufb_otfs_v2','results','r2_2d_confirm_ber');
    addpath(audit,fullfile(pkg,'02_shared_lib'),fullfile(pkg,'04_ber_confirmation'), ...
        fullfile(pkg,'03_axis_specific_pufb'),fullfile(pkg,'01_channel_generator'));
    if ~isfolder(output),mkdir(output);end
    if strcmp(stage,'PREPARE'),result=prepare(root,audit,output,source);return;end
    verify(output);p=jsondecode(fileread(fullfile(output,'r22e_preregistration.json')));p.eta_grid=p.eta_grid(:).';p.development_snr=p.development_snr(:).';p.ber_snr=p.ber_snr(:).';
    physical=jsondecode(fileread(fullfile(root,'pufb_otfs_v2','audit','r2_0a_abid_true2d_pilot','r20a_preregistration.json')));
    if strcmp(stage,'DEVELOP'),result=develop(output,p,physical);return;end
    if strcmp(stage,'FREEZE')
        assert(~isfile(fullfile(output,'r22e_pretest_freeze_receipt.json'))&&isempty(dir(fullfile(output,'test_channels','*.mat'))));
        inputs={fullfile(output,'r22e_supports.mat');fullfile(output,'r22e_support_hashes.csv');fullfile(output,'r22e_eta_selection.csv'); ...
            fullfile(output,'r22e_eta_freeze.mat');fullfile(output,'r22e_val_confirmation.csv');fullfile(output,'r22e_code_hashes.csv'); ...
            fullfile(output,'r22e_train_bank_hashes.csv');fullfile(output,'r22e_support_temporal_receipt.json')};
        writetable(r22d_model('hashes',inputs),fullfile(output,'r22e_pretest_input_hashes.csv'));
        receipt=struct('TRANSFORMS_FROZEN',true,'SUPPORTS_FROZEN',true,'ETAS_FROZEN',true,'K_FROZEN',true,'CODE_FROZEN',true, ...
            'BER_RULES_FROZEN',true,'BOOTSTRAP_RULES_FROZEN',true,'FIGURE_REPORT_CODE_FROZEN',true,'TEST_GENERATED',false,'TEST_ACCESSED',false, ...
            'utc',char(datetime('now','TimeZone','UTC')),'input_hashes_sha256',r22a_io('hash',fullfile(output,'r22e_pretest_input_hashes.csv')));
        r22a_io('json',fullfile(output,'r22e_pretest_freeze_receipt.json'),receipt);
        assert(isequal(receipt,jsondecode(fileread(fullfile(output,'r22e_pretest_freeze_receipt.json')))));
        result=receipt;return
    end
    receipt=jsondecode(fileread(fullfile(output,'r22e_pretest_freeze_receipt.json')));
    assert(receipt.ETAS_FROZEN&&receipt.SUPPORTS_FROZEN&&receipt.CODE_FROZEN&&~receipt.TEST_ACCESSED&&~receipt.TEST_GENERATED);
    assert(strcmp(receipt.input_hashes_sha256,r22a_io('hash',fullfile(output,'r22e_pretest_input_hashes.csv'))));
    r22a_io('verify',fullfile(output,'r22e_pretest_input_hashes.csv'));
    loaded=load(fullfile(output,'r22e_eta_freeze.mat'),'etas','common_eta');etas=loaded.etas;common_eta=loaded.common_eta;
    loaded=load(fullfile(output,'r22e_supports.mat'),'methods');methods=loaded.methods;
    if strcmp(stage,'TEST')
        result=exact_test(methods,etas,common_eta,p,physical,output);
    elseif strcmp(stage,'BER')
        r22e_ber(methods,etas,p,physical,output);result=r22e_finish(methods,p,output);
        status(output,'COMPLETED');artifacts=r22a_io('manifest',output,audit);r22a_io('verify',fullfile(output,'artifact_manifest.csv'));
        fprintf('R22E COMPLETE %s; %d artifacts\n',result.PRIMARY_INTERPRETATION,height(artifacts));
    else
        error('R22E:Stage','Use PREPARE, DEVELOP, FREEZE, TEST or BER.');
    end
    verify(output);r22a_io('verify',fullfile(output,'r22e_pretest_input_hashes.csv'));
end

function p=prepare(root,audit,output,source)
    assert(isempty(dir(fullfile(output,'train_channels','*.mat')))&&~isfile(fullfile(output,'r22e_temporal_freeze_receipt.json')));
    r22a_io('verify',fullfile(source,'artifact_manifest.csv'));r22a_io('verify',fullfile(source,'r22dcb_history_hashes.csv'));
    history=r22e_math('read',fullfile(source,'artifact_manifest.csv'));history=history(:,{'AbsolutePath','SHA256'});
    earlier=r22e_math('read',fullfile(source,'r22dcb_history_hashes.csv'));history=[history;earlier(:,{'AbsolutePath','SHA256'})];
    writetable(history,fullfile(output,'r22e_history_hashes.csv'));
    loaded=load(fullfile(source,'r22dcb_frozen_methods.mat'),'methods');methods=loaded.methods;
    for index=1:4
        checked=r22d_model('transform',methods(index).AlphaNu,methods(index).AlphaTau);assert(isequal(checked.U,methods(index).U));
        methods(index).Mask=[];
    end
    save(fullfile(output,'r22e_frozen_transforms.mat'),'methods','-v7.3');
    p=struct('campaign','R2.2E_SINGLE_NORMALIZED_LOADING_FLOOR','M_FB',2,'K',33178,'Lw',4,'Jnu',2,'Jtau',3, ...
        'train_count',60,'val_count',40,'test_count',200,'development_snr',20:36,'ber_snr',0:36, ...
        'eta_grid',[0,1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1,3e-1,1], ...
        'velocity_mps',[50/3.6,100/3.6,500/3.6,340.3,1701.5], ...
        'velocity_labels',{{'50 km/h','100 km/h','500 km/h','Mach 1','Mach 5'}}, ...
        'seeds',struct('TRAIN',161100000,'VAL',161200000,'TEST',161300000,'DATA',162100000,'NOISE',163100000,'BOOTSTRAP',164100000), ...
        'seed_rule','Physical base+1000*velocity+channel; DATA/NOISE base+100000*velocity+1000*SNR+channel', ...
        'support_rule','Section7: NEW TRAIN_RX-only topK meanabs(U^HGU)^2 masks; exactfixedK and historicalranking; masks frozen before eta search, never updated later. No historicalchannel used.', ...
        'selection_rule','Min equalSNR/equalvelocity TRAIN MSE,within relative1e-4 of minimum choose lowesteta; VAL confirmation only/no retuning', ...
        'common_eta_rule','Min equalmean TRAIN risk over4representations with identicalgrid/tierule', ...
        'receiver_rule','lambda=sigma2+eta*norm(Hs,F)^2/n; onlyHs,sigma2,frozeneta operational. H/E offline objective/diagnostics only', ...
        'val_regions',struct('LOW',20:24,'MID',25:29,'HIGH',30:36),'bootstrap_replicates',10000, ...
        'bootstrap_rule','Paired physicalchannel withinvelocity; equalvelocity averages; noindependentSNR/frame resampling', ...
        'ber_channels_per_velocity',200,'ber_target_errors',500,'ber_cap_bits',20000000,'ber_batch_rounds',16, ...
        'ber_stop_rule','COMMONframes AXIS_REG/SHARED_REG/D4_REG/AXIS_ETA0 all>=500errors OR floor(2e7/(576*200))rounds; R22A_REG sameframes secondary', ...
        'ber_targets',[1e-2,1e-3,1e-4],'ber_modulation','GrayuncodedQPSK as historicalimplementation; paired bits/symbols/physicalnoise all5methods', ...
        'ber_aggregate','Equal5velocity count-ratios; pooledraw separate descriptive', ...
        'crossing_rule','Firstdescending adjacentpositive measuredBERbracket on0:36; require>=500original aggregateerrors at bothpoints, noextrapolation/fitting', ...
        'gain_support_rule','Atleast95% pairedphysicalchannelbootstrapcrossingsresolved AND percentileCI lower>0; conditionalCI and resolutioncounts explicit', ...
        'shape_rule','OLS log10positive measuredBER at30:36; >=3positive pointsrequired; <=-.03 CLEAR_DESCENT;(-.03,.03) APPROX_FLAT;>=.03 UPTURN; elseUNRESOLVED', ...
        'numerical_rule','First10TESTchannelspervelocity,all10exactvariants at20/30/36: compareSVD risk toCholeskyresidual<=1e-10*(1+MSE),backward<=1e-12; fullMSEmonotonic tolerance1e-10', ...
        'mechanism_rule','AxisREG-minusETA0 TESTmean MSE upperCI<0 AND HIGH30:36 WSE_F,T_OMIT upperCI<0 AND numericalPASS', ...
        'upturn_rule','REG slope<ETA0slope; pairedHIGH30:36BERdiff upperCI<0; MSEbenefit/WSEreduction/numericalPASS. Strong interpretation requiresREG CLEAR_DESCENT', ...
        'robust_representation_rule','AxisREG-minusSHAREDREG MSE upperCI<0 AND AxisCOMMON-minusSHAREDCOMMON upperCI<0; bothdirectsamereceiverfamily', ...
        'interpretation_precedence','Etaaxis0 NO_SUPPORTED_NONZERO; TESTMSEfailure REGULARIZATION_EFFECT_NOT_CONFIRMED_OUT_OF_SAMPLE; mechanismandupturncorrectedandCLEAR_DESCENT MISMATCH_AWARE_DIAGONAL_LOADING_CORRECTS_HIGH_SNR_SPARSE_RECEIVER_DEGRADATION; supportedBERwithoutmechanism BER_IMPROVEMENT_WITHOUT_CONFIRMED_MISMATCH_MECHANISM; otherwise SECOND_ORDER_RECEIVER_GAIN_WITHOUT_RESOLVED_BER_GAIN', ...
        'BER_USED_FOR_ETA_SELECTION',false,'TEST_USED_FOR_ETA_SELECTION',false,'HISTORICAL_BER_EXT_USED_FOR_SELECTION',false, ...
        'NEW_TRANSFORM_TRAINED',false,'ALPHA_CHANGED',false,'FILTER_LENGTH_ESCALATED',false,'MANUSCRIPT_MODIFIED',false);
    path=fullfile(output,'r22e_preregistration.json');
    if isfile(path),assert(isequaln(jsondecode(fileread(path)),jsondecode(jsonencode(p))));else,r22a_io('json',path,p);end
    grid=table((1:14)',p.eta_grid(:),'VariableNames',{'CandidateID','Eta'});writetable(grid,fullfile(output,'r22e_eta_grid.csv'));
    manifest=table();
    for split={'train','val','test'}
        name=split{1};count=p.([name,'_count']);v=repelem((1:5)',count);c=repmat((1:count)',5,1);seed=p.seeds.(upper(name))+1000*v+c;
        T=table(repmat(string(upper(name)),5*count,1),v,c,seed,'VariableNames',{'Split','VelocityIndex','ChannelID','Seed'});manifest=[manifest;T]; %#ok<AGROW>
    end
    old_files=dir(fullfile(root,'pufb_otfs_v2','results','**','*manifest.csv'));
    for index=1:numel(old_files)
        path=fullfile(old_files(index).folder,old_files(index).name);if startsWith(path,[output,filesep]),continue;end
        T=r22e_math('read',path);if ismember('Seed',T.Properties.VariableNames),assert(isempty(intersect(T.Seed,manifest.Seed)));end
    end
    writetable(manifest,fullfile(output,'r22e_seed_manifest.csv'));writetable(manifest,fullfile(output,'r22e_split_manifest.csv'));
    r22a_io('json',fullfile(output,'r22e_kernel_preflight.json'),r22e_math('selftest'));r22e_finish('SELFTEST',p,output);
    files=dir(fullfile(audit,'*.m'));paths=cell(numel(files),1);
    for index=1:numel(files)
        paths{index}=fullfile(audit,files(index).name);d=checkcode(paths{index},'-id');
        for item=1:numel(d),fprintf('%s:%d %s\n',files(index).name,d(item).line,d(item).message);end
        assert(isempty(d),'R22E:Checkcode','Pre-generation diagnostic.');
    end
    paths=[paths;{which('r22d_model');which('r22dcb_math');which('r22a_io');which('fir_mallat_r22a');which('generate_channel_r20a'); ...
        fullfile(root,'pufb_otfs_v2','audit','r2_0a_abid_true2d_pilot','r20a_preregistration.json')}];
    writetable(r22d_model('hashes',paths),fullfile(output,'r22e_code_hashes.csv'));
    inputs={fullfile(source,'r22dcb_frozen_methods.mat');fullfile(output,'r22e_frozen_transforms.mat');fullfile(output,'r22e_preregistration.json'); ...
        fullfile(output,'r22e_seed_manifest.csv');fullfile(output,'r22e_split_manifest.csv');fullfile(output,'r22e_eta_grid.csv');fullfile(output,'r22e_history_hashes.csv')};
    writetable(r22d_model('hashes',inputs),fullfile(output,'r22e_input_hashes.csv'));
    receipt=struct('SPLITS_FROZEN',true,'SEEDS_FROZEN',true,'ETA_GRID_FROZEN',true,'TRANSFORMS_FROZEN',true,'CODE_FROZEN',true, ...
        'NO_CHANNELS_GENERATED',true,'TEST_INACCESSIBLE',true,'utc',char(datetime('now','TimeZone','UTC')), ...
        'input_hashes_sha256',r22a_io('hash',fullfile(output,'r22e_input_hashes.csv')),'code_hashes_sha256',r22a_io('hash',fullfile(output,'r22e_code_hashes.csv')));
    r22a_io('json',fullfile(output,'r22e_temporal_freeze_receipt.json'),receipt);verify(output);status(output,'PREREGISTERED_BEFORE_CHANNELS');
end

function result=develop(output,p,physical)
    assert(~isfile(fullfile(output,'r22e_pretest_freeze_receipt.json')),'R22E:Frozen','No selection afterTESTfreeze.');
    train=bank('TRAIN',p,physical,output);loaded=load(fullfile(output,'r22e_frozen_transforms.mat'),'methods');methods=loaded.methods;
    support_path=fullfile(output,'r22e_supports.mat');
    if isfile(support_path),loaded=load(support_path,'methods');methods=loaded.methods;
    else
        rows=cell(4,3);
        for method=1:4
            methods(method).Mask=r22d_model('support',train,methods(method).U,p.K);rows(method,:)={methods(method).Name,p.K,r22dcb_math('arrayhash',uint8(methods(method).Mask(:)))};
        end
        save(support_path,'methods','-v7.3');writetable(cell2table(rows,'VariableNames',{'Method','RetainedEntries','MaskSHA256'}),fullfile(output,'r22e_support_hashes.csv'));
        files=dir(fullfile(output,'train_channels','*.mat'));paths=arrayfun(@(f)fullfile(f.folder,f.name),files,'UniformOutput',false);
        writetable(r22d_model('hashes',paths),fullfile(output,'r22e_train_bank_hashes.csv'));
        r22a_io('json',fullfile(output,'r22e_support_temporal_receipt.json'),struct('SUPPORTS_FROZEN_BEFORE_ETA_SEARCH',true, ...
            'TRAIN_RX_ONLY',true,'VAL_OPENED',false,'TEST_OPENED',false,'K',p.K,'SupportMAT_SHA256',r22a_io('hash',support_path)));
    end
    train_risk=development_risk(train,methods,p.eta_grid,p,output,'train');risk=squeeze(mean(train_risk,[1,2]));
    rows=cell(56,3);cursor=0;etas=zeros(1,4);
    for method=1:4
        index=r22e_math('select',risk(:,method),p.eta_grid);etas(method)=p.eta_grid(index);
        for candidate=1:14,cursor=cursor+1;rows(cursor,:)={methods(method).Name,p.eta_grid(candidate),risk(candidate,method)};end
    end
    common_index=r22e_math('select',mean(risk,2),p.eta_grid);common_eta=p.eta_grid(common_index);
    writetable(cell2table(rows,'VariableNames',{'Method','Eta','TrainJrx'}),fullfile(output,'r22e_train_eta_objective.csv'));
    selection=table([string({methods.Name}).';"COMMON"],[etas(:);common_eta], 'VariableNames',{'Method','EtaStar'});
    writetable(selection,fullfile(output,'r22e_eta_selection.csv'));save(fullfile(output,'r22e_eta_freeze.mat'),'etas','common_eta','-v7.3');
    val=bank('VAL',p,physical,output);val_values=zeros(200,17,3,4);
    for method=1:4
        current=development_risk(val,methods(method),[0,etas(method),common_eta],p,output,sprintf('val_m%d',method));val_values(:,:,:,method)=current;
    end
    rows=cell(4*4,7);cursor=0;regions={20:36,20:24,25:29,30:36};names={'AGGREGATE','LOW','MID','HIGH'};
    for method=1:4
        for region=1:4
            points=ismember(p.development_snr,regions{region});left=mean(val_values(:,points,1,method),2);right=mean(val_values(:,points,2,method),2);
            stats=r22e_math('bootstrap',right-left,40,p.seeds.BOOTSTRAP);cursor=cursor+1;
            rows(cursor,:)={methods(method).Name,string(names{region}),mean(left),mean(right),stats.Mean,stats.Lower,stats.Upper};
        end
    end
    validation=cell2table(rows,'VariableNames',{'Method','Region','Eta0MSE','EtaStarMSE','MeanDifference','CI95Lower','CI95Upper'});
    writetable(validation,fullfile(output,'r22e_val_confirmation.csv'));save(fullfile(output,'r22e_development_context.mat'),'risk','validation','-v7.3');
    status(output,'ETAS_TRAIN_SELECTED_AND_VAL_CONFIRMED_NO_RETUNING');result=selection;
end

function risk=development_risk(values,methods,etas,p,output,label)
    directory=fullfile(output,[label,'_risk']);if ~isfolder(directory),mkdir(directory);end
    risk=zeros(size(values,3),17,numel(etas),numel(methods));
    for channel=1:size(values,3)
        path=fullfile(directory,sprintf('channel_%04d.mat',channel));
        if isfile(path),loaded=load(path,'measured');measured=loaded.measured;
        else
            measured=zeros(17,numel(etas),numel(methods));
            for method=1:numel(methods)
                H=methods(method).U'*values(:,:,channel)*methods(method).U;geometry=r22e_math('geometry',H,methods(method).Mask);
                measured(:,:,method)=r22e_math('risk',geometry,etas,p.development_snr,false).';
            end
            save(path,'measured');
        end
        risk(channel,:,:,:)=measured;
        if mod(channel,50)==0,fprintf('R22E %s exact risk %d/%d\n',label,channel,size(values,3));end
    end
end

function variants=make_variants(methods,etas,common_eta)
    indices=[3,3,2,2,1,4,1,2,3,4];names=["AXIS_ETA0","AXIS_REG","SHARED_ETA0","SHARED_REG","D4_REG","R22A_REG", ...
        "D4_COMMON","SHARED_COMMON","AXIS_COMMON","R22A_COMMON"];
    parameters=[0,etas(3),0,etas(2),etas(1),etas(4),repmat(common_eta,1,4)];
    variants=struct('Name',cell(1,10),'MethodIndex',cell(1,10),'Eta',cell(1,10));
    for index=1:10,variants(index).Name=names(index);variants(index).MethodIndex=indices(index);variants(index).Eta=parameters(index);end
    assert(numel(methods)==4);
end

function result=exact_test(methods,etas,common_eta,p,physical,output)
    variants=make_variants(methods,etas,common_eta);save(fullfile(output,'r22e_test_variants.mat'),'variants');
    values=bank('TEST',p,physical,output);rows=cell(170000,7);diagnostic_rows=cell(34000,15);cursor=0;diagnostic_cursor=0;
    full_rows=cell(68000,5);full_cursor=0;maximum_error=0;maximum_backward=0;monotonic=true;
    for channel=1:1000
        path=fullfile(output,sprintf('test_metric_%04d.mat',channel));
        if isfile(path),loaded=load(path,'record');record=loaded.record;
        else
            record=struct('MSE',zeros(10,17),'FullMSE',zeros(4,17),'AxisDiagnostics',zeros(2,17,11),'DirectError',0,'Backward',0);
            for method=1:4
                H=methods(method).U'*values(:,:,channel)*methods(method).U;geometry=r22e_math('geometry',H,methods(method).Mask);
                indices=find([variants.MethodIndex]==method);parameters=[variants(indices).Eta];
                measured=r22e_math('risk',geometry,parameters,p.development_snr,false);record.MSE(indices,:)=measured;
                singular=svd(H);record.FullMSE(method,:)=mean(10.^(-p.development_snr/10)./(singular.^2+10.^(-p.development_snr/10)),1);
                if method==3,record.AxisDiagnostics=r22e_math('risk',geometry,[0,etas(3)],p.development_snr,true);end
                if mod(channel-1,200)<10
                    for index=indices
                        for point=[1,11,17]
                            variance=10^(-p.development_snr(point)/10);W=r22e_math('receiver',geometry.Hs,variance,variants(index).Eta);
                            direct=(norm(W*H-eye(288),'fro')^2+variance*norm(W,'fro')^2)/288;
                            record.DirectError=max(record.DirectError,abs(direct-record.MSE(index,point)));
                            loading=variance+variants(index).Eta*geometry.Q;A=geometry.Hs'*geometry.Hs+loading*eye(288);
                            residual=norm(A*W-geometry.Hs','fro')/((geometry.Singular(1)^2+loading)*norm(W,'fro')+norm(geometry.Hs,'fro'));
                            record.Backward=max(record.Backward,residual);
                        end
                    end
                end
            end
            assert(max(abs(record.FullMSE-record.FullMSE(1,:)),[],'all')<1e-12&&record.DirectError<1e-10&&record.Backward<1e-12);
            save(path,'record','-v7.3');
        end
        maximum_error=max(maximum_error,record.DirectError);maximum_backward=max(maximum_backward,record.Backward);
        monotonic=monotonic&&all(diff(record.FullMSE,1,2)<=1e-10,'all');velocity=ceil(channel/200);local=mod(channel-1,200)+1;
        for index=1:10
            for point=1:17,cursor=cursor+1;rows(cursor,:)={variants(index).Name,velocity,local,p.development_snr(point),variants(index).Eta,p.K,record.MSE(index,point)};end
        end
        for index=1:2
            name="AXIS_ETA0";eta=0;if index==2,name="AXIS_REG";eta=etas(3);end
            for point=1:17
                diagnostic_cursor=diagnostic_cursor+1;diagnostic_rows(diagnostic_cursor,:)=[{name,velocity,local,p.development_snr(point)},num2cell(reshape(record.AxisDiagnostics(index,point,:),1,11))];
            end
            assert(eta>=0);
        end
        for method=1:4
            for point=1:17,full_cursor=full_cursor+1;full_rows(full_cursor,:)={methods(method).Name,velocity,local,p.development_snr(point),record.FullMSE(method,point)};end
        end
        if mod(channel,50)==0,fprintf('R22E TEST exact metrics %d/1000\n',channel);end
    end
    mse=cell2table(rows,'VariableNames',{'Variant','VelocityIndex','ChannelID','SNRdB','Eta','RetainedEntries','MSE'});
    diagnostic=cell2table(diagnostic_rows,'VariableNames',{'Variant','VelocityIndex','ChannelID','SNRdB','MSE','T_REG','T_OMIT','T_CROSS','T_NOISE', ...
        'W_F','WSE_F','CondA','LambdaEff','Qs','MismatchPenalty'});
    writetable(mse,fullfile(output,'r22e_test_exact_mse.csv'));writetable(diagnostic,fullfile(output,'r22e_test_mismatch_decomposition.csv'));
    writetable(diagnostic(:,{'Variant','VelocityIndex','ChannelID','SNRdB','W_F','WSE_F','CondA','LambdaEff','Qs'}),fullfile(output,'r22e_test_receiver_gain.csv'));
    writetable(mse(endsWith(mse.Variant,"COMMON"),:),fullfile(output,'r22e_common_eta_control.csv'));
    writetable(cell2table(full_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','FullMSE'}),fullfile(output,'r22e_full_mse.csv'));
    result=struct('NUMERICAL_SOLVER_PASS',maximum_error<1e-10&&maximum_backward<1e-12,'FULL_MSE_MONOTONIC_PASS',monotonic, ...
        'MaxDirectSolveMSEError',maximum_error,'MaxBackwardResidual',maximum_backward);
    r22a_io('json',fullfile(output,'r22e_numerical_controls.json'),result);status(output,'EXACT_TEST_COMPLETE_ETAS_UNCHANGED');
end

function values=bank(split,p,physical,output)
    manifest=r22e_math('read',fullfile(output,'r22e_split_manifest.csv'));manifest=manifest(manifest.Split==string(split),:);
    values=complex(zeros(288,288,height(manifest)));directory=fullfile(output,[lower(split),'_channels']);if ~isfolder(directory),mkdir(directory);end
    for index=1:height(manifest)
        path=fullfile(directory,sprintf('channel_%04d.mat',index));
        if isfile(path),loaded=load(path,'G','seed');G=loaded.G;assert(loaded.seed==manifest.Seed(index));
        else,seed=manifest.Seed(index);G=generate_channel_r20a(24,12,p.velocity_mps(manifest.VelocityIndex(index)),seed,physical);save(path,'G','seed','-v7.3');end
        values(:,:,index)=G;
    end
end

function verify(output)
    for name={'r22e_code_hashes.csv','r22e_input_hashes.csv','r22e_history_hashes.csv'},r22a_io('verify',fullfile(output,name{1}));end
end

function status(output,label)
    r22a_io('json',fullfile(output,'r22e_execution_status.json'),struct('status',label,'NEW_TRANSFORM_TRAINED',false, ...
        'TEST_USED_FOR_ETA_SELECTION',false,'BER_USED_FOR_ETA_SELECTION',false,'HISTORICAL_BER_EXT_USED_FOR_SELECTION',false,'MANUSCRIPT_MODIFIED',false));
end