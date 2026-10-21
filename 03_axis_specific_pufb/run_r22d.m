function result=run_r22d(root,stage)
    pkg=fileparts(mfilename('fullpath'));pkg=fileparts(pkg);
    audit=fullfile(pkg,'03_axis_specific_pufb');
    output=fullfile(root,'pufb_otfs_v2','results','r2_2d_axis_specific');
    helper=fullfile(pkg,'02_shared_lib');
    physical_code=fullfile(pkg,'01_channel_generator');
    physical_cfg=fullfile(root,'pufb_otfs_v2','audit','r2_0a_abid_true2d_pilot');
    addpath(audit,helper,physical_code,fullfile(root,'pufb_otfs_v2','audit','r2_2c_confirm_strict'));
    if ~isfolder(output),mkdir(output);end
    if strcmp(stage,'PREPARE')
        result=prepare(root,audit,output,helper,physical_code);return
    end
    p=jsondecode(fileread(fullfile(output,'r22d_preregistration.json')));
    physical=jsondecode(fileread(fullfile(physical_cfg,'r20a_preregistration.json')));
    r22a_io('verify',fullfile(output,'r22d_history_preservation_receipt.csv'));
    r22a_io('verify',fullfile(output,'r22d_code_hashes.csv'));
    if strcmp(stage,'TRAIN')
        result=train_select(root,output,p,physical);return
    end
    if strcmp(stage,'FREEZE')
        result=freeze_test(audit,output,p);return
    end
    assert(strcmp(stage,'TEST'),'R22D:Stage','Use PREPARE, TRAIN, FREEZE or TEST.');
    receipt=jsondecode(fileread(fullfile(output,'r22d_temporal_freeze_receipt.json')));
    assert(receipt.CANDIDATES_FROZEN&&receipt.SUPPORTS_FROZEN&&receipt.COMPARATORS_FROZEN&&receipt.CODE_FROZEN ...
        &&~receipt.TEST_GENERATED&&~receipt.TEST_ACCESSED,'R22D:TemporalGate','Invalid temporal gate.');
    assert(strcmp(receipt.input_hashes_sha256,r22a_io('hash',fullfile(output,'r22d_test_input_hashes.csv'))));
    r22a_io('verify',fullfile(output,'r22d_test_input_hashes.csv'));
    loaded=load(fullfile(output,'r22d_supports.mat'),'methods');methods=loaded.methods;
    write_status(output,'TEST_OPENED_AFTER_COMPLETE_FREEZE',true,true);
    [mse,frob,sparsity,full,stability]=r22d_test(methods,p,physical,output);
    result=r22d_finish(methods,p,mse,frob,sparsity,full,stability,output);
    r22a_io('verify',fullfile(output,'r22d_code_hashes.csv'));
    r22a_io('verify',fullfile(output,'r22d_test_input_hashes.csv'));
    r22a_io('verify',fullfile(output,'r22d_history_preservation_receipt.csv'));
    write_status(output,'COMPLETED',true,true);
    manifest=r22a_io('manifest',output,audit);r22a_io('verify',fullfile(output,'artifact_manifest.csv'));
    fprintf('R22D COMPLETE %s; %d artifacts\n',result.PRIMARY_INTERPRETATION,height(manifest));
end

function p=prepare(root,audit,output,helper,physical_code)
    prereg_path=fullfile(output,'r22d_preregistration.json');
    physical_cfg=fullfile(root,'pufb_otfs_v2','audit','r2_0a_abid_true2d_pilot');
    assert(isempty(dir(fullfile(output,'train_channels','*.mat'))),'R22D:Prepared','No re-preparation after opening TRAIN.');
    canonical=jsondecode(fileread(fullfile(root,'pufb_otfs_v2','audit','r2_2c_receiver_aware_fixed_budget','r22c_preregistration.json')));
    assert(startsWith(canonical.alpha_canonical_interval,'[0,pi)'));
    history_dirs={'r2_0a_abid_true2d_pilot','r2_2a_abid_fir_sparse','r2_2b_receiver_useful_sparsity_audit', ...
        'r2_2c_receiver_aware_fixed_budget','r2_2c_confirm_strict'};history=table();
    for index=1:numel(history_dirs)
        file=fullfile(root,'pufb_otfs_v2','results',history_dirs{index},'artifact_manifest.csv');
        r22a_io('verify',file);T=r22d_model('read',file);history=[history;T(:,{'AbsolutePath','SHA256'})]; %#ok<AGROW>
    end
    writetable(history,fullfile(output,'r22d_history_preservation_receipt.csv'));
    old_star=load(fullfile(root,'pufb_otfs_v2','results','r2_2a_abid_fir_sparse','r22a_pufb2d_star.mat'),'alpha');
    assert(abs(old_star.alpha-.7967306572)<5e-11);
    p=struct('campaign','R2.2D_AXIS_SPECIFIC_RECEIVER_AWARE','physical_scenario','ABID_EVA', ...
        'M',24,'N',12,'L',288,'modulation_order',4,'carrier_frequency_hz',4e12,'subcarrier_spacing_hz',15e6, ...
        'filter_length',4,'J_nu',2,'J_tau',3,'K',33178,'ALPHA_NU_DOMAIN','[0,pi)','ALPHA_TAU_DOMAIN','[0,pi)', ...
        'screen_per_velocity',20,'train_per_velocity',40,'val_per_velocity',20,'test_per_velocity',200, ...
        'velocity_mps',[50/3.6,100/3.6,500/3.6,340.3,1701.5], ...
        'velocity_labels',{{'50 km/h','100 km/h','500 km/h','Mach 1','Mach 5'}}, ...
        'seeds',struct('TRAIN',121100000,'VAL',121200000,'TEST',121300000,'BOOTSTRAP',121400000,'SYNTHETIC',119001), ...
        'seed_rule','base+1000*velocity_index+local_channel_id', ...
        'screen_channels','Local channel IDs1:20 in each TRAIN velocity; only these generated/opened for coarse shared search.', ...
        'screen_snr_db',[6,15,24],'snr_db',3:3:24,'strata',struct('LOW',[3,6,9],'MID',[12,15,18],'HIGH',[21,24]), ...
        'shared_coarse_points',61,'shared_neighborhoods',3,'shared_refine_points',21, ...
        'shared_neighborhood_rule','Periodic one-dimensional local minima ranked by screen J_RX, then fill from ranked samples; centers >=3*pi/61 apart circularly; refinement radius pi/61.', ...
        'axis_grid_points',17,'axis_neighborhoods',4,'axis_refine_points',9, ...
        'axis_neighborhood_rule','Toroidal 8-neighbor minima on17x17 plus inserted points, ranked by screen J_RX, then fill ranked samples. Accept center only if toroidal separation >=2*pi/17 on at least one axis from every accepted center.', ...
        'axis_refine_rule','9x9 centered at each selected coarse/inserted pair, radius pi/(2*17) independently per axis; canonical modulo pi; round14decimal for duplicate removal only.', ...
        'final_axis_mesh',pi/(17*8),'parameter_separation_multiple',2, ...
        'separation_rule','Absolute reported alpha_tau-alpha_nu must exceed2*final_axis_mesh; also report circular separation as a caveat near canonical boundaries.', ...
        'candidate_ties','J_RX then alpha_nu then alpha_tau ascending; no TEST tie break', ...
        'frozen_shared_candidates',3,'frozen_axis_candidates',5,'validation_access_count',1, ...
        'objective','Mean of LOW/MID/HIGH excess-MSE ratios to D4 baseline on identical stage TRAIN bank, support budget and SNR grid; equal velocity weighting.', ...
        'val_normalization','Same fixed full-TRAIN D4 stratum denominators for both arms; never fit VAL baseline or support.', ...
        'primary_test_weighting','Per channel mean LOW + mean MID + mean HIGH divided by3; then equal velocity means. Not uniform over8SNR points.', ...
        'bootstrap_replicates',10000,'bootstrap_unit','Paired physical channels; resample200 channels within each velocity independently; equally average five velocity means; never resample SNR independently.', ...
        'ci_rule','Percentile sorted order statistics ceil(.025*10000),ceil(.975*10000); SUPPORTED upper<0; NOT_SUPPORTED lower>=0; otherwise INCONCLUSIVE', ...
        'stress_velocity',5,'stress_snr_db',[21,24],'structural_tolerance',1e-11,'full_invariance_tolerance',1e-12, ...
        'support_rule','Top33178 mean abs(U^H*G*U)^2 entries on scoring TRAIN ensemble; support fixed across its channels. Freeze full-TRAIN masks before VAL.', ...
        'stability_diagnostic_rule','Oracle topK diagnostic only after TEST; first40 channels per velocity for all within-velocity pairs and cross-velocity pairs against each other velocity. Never used in receivers or selection.', ...
        'strong_regime_dependence_rule','Aggregate supported AND at least one velocity or SNR-stratum lower CI>=0 (significant adverse difference); secondary intervals unadjusted.', ...
        'interpretation_precedence','F significant regime dependence; A strong specialization; C supported but separation unresolved; B supported without assignment support; G null and circular distance<=mesh; E null same-bank but historical R22C contrast supported; D remaining null.', ...
        'historical_r22c_alpha',.78797323934301,'historical_r22a_alpha',old_star.alpha, ...
        'BER_RUN',false,'PCG_RUN',false,'K99_USED_FOR_SELECTION',false,'BER_USED_FOR_TRAINING',false, ...
        'GRAM_PENALTY_USED',false,'FILTER_LENGTH_ESCALATED',false,'NONSEPARABLE_TRANSFORM_USED',false,'MANUSCRIPT_MODIFIED',false, ...
        'metric_implementation','Exact eigenbasis solution of regularized Gram; no eigenmode truncation, no proxy. Validated against Cholesky residual at preflight. TEST uses direct Cholesky residual on sparse and full channels.', ...
        'doppler','Unmodified generate_channel_r20a; k=nu*N*T');
    prior_files=dir(fullfile(root,'pufb_otfs_v2','results','**','*manifest.csv'));prior_seeds=[];
    for index=1:numel(prior_files)
        path=fullfile(prior_files(index).folder,prior_files(index).name);
        if startsWith(path,[output,filesep]),continue;end
        T=r22d_model('read',path);
        if ismember('Seed',T.Properties.VariableNames),prior_seeds=[prior_seeds;T.Seed];end %#ok<AGROW>
    end
    for split={'train','val','test'}
        label=split{1};count=p.([label,'_per_velocity']);velocity=repelem((1:5)',count);channel=repmat((1:count)',5,1);
        seed=p.seeds.(upper(label))+1000*velocity+channel;assert(isempty(intersect(seed,prior_seeds)));
        T=table(velocity,channel,repelem(p.velocity_mps(:),count),seed,false(5*count,1), ...
            'VariableNames',{'VelocityIndex','ChannelID','VelocityMps','Seed','ChannelGenerated'});
        writetable(T,fullfile(output,['r22d_',label,'_manifest.csv']));
    end
    if isfile(prereg_path)
        previous=jsondecode(fileread(prereg_path));
        assert(isequaln(previous,jsondecode(jsonencode(p))),'R22D:PreregChanged','Existing preregistration differs.');
    else
        r22a_io('json',prereg_path,p);
    end
    metric_check=r22d_model('selftest');
    r22a_io('json',fullfile(output,'r22d_preflight_metric.json'),metric_check);
    r22d_finish('SELFTEST',p,[],[],[],[],[],output);
    files=dir(fullfile(audit,'*.m'));paths=cell(numel(files),1);
    for index=1:numel(files)
        path=fullfile(audit,files(index).name);paths{index}=path;d=checkcode(path,'-id');
        for diagnostic=1:numel(d),fprintf('%s:%d %s\n',files(index).name,d(diagnostic).line,d(diagnostic).message);end
        assert(isempty(d),'R22D:Checkcode','Source diagnostics before channels.');
    end
    paths=[paths;{fullfile(helper,'fir_mallat_r22a.m');fullfile(helper,'r22a_io.m'); ...
        fullfile(root,'pufb_otfs_v2','audit','r2_2c_confirm_strict','confirm_exact_mse.m'); ...
        fullfile(physical_code,'generate_channel_r20a.m');fullfile(physical_cfg,'r20a_preregistration.json');prereg_path}];
    writetable(r22d_model('hashes',paths),fullfile(output,'r22d_code_hashes.csv'));
    write_status(output,'PREREGISTERED_BEFORE_NEW_BANKS',false,false);
end

function result=train_select(root,output,p,physical)
    assert(~isfile(fullfile(output,'r22d_temporal_freeze_receipt.json')),'R22D:AlreadyFrozen','No training after TEST freeze.');
    screen=bank('train',p.screen_per_velocity,p,physical,output);
    screen_singular=singular_bank(screen);
    d4_screen=r22d_model('score',[pi/3,pi/3],screen,screen_singular,[],p.screen_snr_db,p,[]);
    alpha=(0:60)'*pi/61;
    coarse=evaluate_pairs([alpha,alpha],screen,screen_singular,d4_screen.Raw,p.screen_snr_db,p,output,'shared_coarse');
    writetable(coarse,fullfile(output,'r22d_shared_coarse_search.csv'));
    shared_centers=shared_neighborhoods(coarse,p);shared_pairs=[];
    for center=shared_centers(:).'
        values=mod(center+linspace(-pi/61,pi/61,21)',pi);shared_pairs=[shared_pairs;values,values]; %#ok<AGROW>
    end
    shared_pairs=unique(round(shared_pairs,14),'rows');
    train=bank('train',p.train_per_velocity,p,physical,output);train_singular=singular_bank(train);
    d4_full=r22d_model('score',[pi/3,pi/3],train,train_singular,[],p.snr_db,p,[]);
    baseline=d4_full.Raw;save(fullfile(output,'r22d_train_baseline.mat'),'baseline');
    shared=evaluate_pairs(shared_pairs,train,train_singular,baseline,p.snr_db,p,output,'shared_refine');
    writetable(shared,fullfile(output,'r22d_shared_refinement.csv'));ranked_shared=sortrows(shared,{'Jrx','AlphaNu','AlphaTau'});
    axis_values=(0:16)'*pi/17;[nu,tau]=ndgrid(axis_values,axis_values);grid_pairs=[nu(:),tau(:)];
    extras=[p.historical_r22c_alpha,p.historical_r22c_alpha;p.historical_r22a_alpha,p.historical_r22a_alpha; ...
        pi/3,pi/3;ranked_shared.AlphaNu(1),ranked_shared.AlphaTau(1)];
    axis_pairs=unique([grid_pairs;extras],'rows','stable');
    global_surface=evaluate_pairs(axis_pairs,screen,screen_singular,d4_screen.Raw,p.screen_snr_db,p,output,'axis_global');
    writetable(global_surface,fullfile(output,'r22d_axis_global_surface.csv'));
    axis_centers=axis_neighborhoods(global_surface,p);refine_pairs=[];
    for center=1:4
        nu_values=mod(axis_centers(center,1)+linspace(-pi/34,pi/34,9),pi);
        tau_values=mod(axis_centers(center,2)+linspace(-pi/34,pi/34,9),pi);
        [nu,tau]=ndgrid(nu_values,tau_values);refine_pairs=[refine_pairs;nu(:),tau(:)]; %#ok<AGROW>
    end
    refine_pairs=unique(round(refine_pairs,14),'rows');
    refined=evaluate_pairs(refine_pairs,train,train_singular,baseline,p.snr_db,p,output,'axis_refine');
    writetable(refined,fullfile(output,'r22d_axis_refinement.csv'));
    surface=[decorate(coarse,"SHARED_SCREEN");decorate(shared,"SHARED_FULL_TRAIN"); ...
        decorate(global_surface,"AXIS_SCREEN");decorate(refined,"AXIS_FULL_TRAIN")];
    writetable(surface,fullfile(output,'r22d_train_surface.csv'));
    ranked_axis=sortrows(refined,{'Jrx','AlphaNu','AlphaTau'});
    shared_selected=ranked_shared(1:3,:);axis_selected=ranked_axis(1:5,:);
    candidate_path=fullfile(output,'r22d_train_candidates_frozen.mat');
    if ~isfile(candidate_path)
        candidates=cell(8,1);
        for index=1:8
            if index<=3,T=shared_selected(index,:);prefix='shared_refine';else,T=axis_selected(index-3,:);prefix='axis_refine';end
            stored=load(candidate_file(output,prefix,T.CandidateID),'candidate');candidates{index}=stored.candidate;
        end
        save(candidate_path,'candidates','shared_selected','axis_selected','shared_centers','axis_centers','-v7.3');
        r22a_io('json',fullfile(output,'r22d_train_candidate_freeze.json'),struct('sha256',r22a_io('hash',candidate_path), ...
            'shared_candidates',3,'axis_candidates',5,'VAL_ACCESSED',false));
    end
    loaded=load(candidate_path,'candidates');candidates=loaded.candidates;
    selection_path=fullfile(output,'r22d_val_selection.csv');
    if ~isfile(selection_path)
        write_status(output,'ALL_TRAIN_CANDIDATES_FROZEN_VAL_OPENED_ONCE',false,false);
        val=bank('val',p.val_per_velocity,p,physical,output);val_singular=singular_bank(val);rows=cell(8,7);
        for index=1:8
            candidate=candidates{index};scored=r22d_model('score',candidate.Pair,val,val_singular,baseline,p.snr_db,p,candidate.Method.Mask);
            arm="AXIS";if index<=3,arm="SHARED";end
            rows(index,:)={index,arm,candidate.Pair(1),candidate.Pair(2),candidate.Jrx,scored.Jrx,scored.Ratios};
        end
        selection=cell2table(rows,'VariableNames',{'CandidateID','Arm','AlphaNu','AlphaTau','TrainJrx','ValJrx','ValStrata'});
        writetable(selection,selection_path);
    else
        selection=r22d_model('read',selection_path);
    end
    shared_val=sortrows(selection(selection.Arm=="SHARED",:),{'ValJrx','AlphaNu','AlphaTau'});
    axis_val=sortrows(selection(selection.Arm=="AXIS",:),{'ValJrx','AlphaNu','AlphaTau'});
    selected_ids=[shared_val.CandidateID(1),axis_val.CandidateID(1)];
    shared_method=candidates{selected_ids(1)}.Method;shared_method.Name="R22D_SHARED_CONTROL";
    axis_method=candidates{selected_ids(2)}.Method;axis_method.Name="R22D_AXIS_SPECIFIC_STAR";
    shared_val_jrx=shared_val.ValJrx(1);axis_val_jrx=axis_val.ValJrx(1);
    save(fullfile(output,'r22d_primary_arms.mat'),'shared_method','axis_method','shared_val_jrx','axis_val_jrx','-v7.3');
    write_arm_freeze(output,'r22d_shared_freeze.json',shared_method,shared_val_jrx);
    write_arm_freeze(output,'r22d_axis_specific_freeze.json',axis_method,axis_val_jrx);
    methods=repmat(shared_method,1,6);methods(4)=shared_method;methods(5)=axis_method;
    pairs=[pi/3,pi/3;p.historical_r22a_alpha,p.historical_r22a_alpha;p.historical_r22c_alpha,p.historical_r22c_alpha; ...
        shared_method.AlphaNu,shared_method.AlphaTau;axis_method.AlphaNu,axis_method.AlphaTau;axis_method.AlphaTau,axis_method.AlphaNu];
    names=["D4_2D_PARTIAL_MALLAT","R2.2A_PUFB2D_STAR","R2.2C_CONFIRMED_STAR","R22D_SHARED_CONTROL","R22D_AXIS_SPECIFIC_STAR","R22D_AXIS_SWAPPED_CONTROL"];
    for index=[1,2,3,6]
        methods(index)=r22d_model('transform',pairs(index,1),pairs(index,2));
        methods(index).Mask=r22d_model('support',train,methods(index).U,p.K);
    end
    for index=1:6,methods(index).Name=names(index);assert(nnz(methods(index).Mask)==33178);end
    old=load(fullfile(root,'pufb_otfs_v2','results','r2_2a_abid_fir_sparse','r22a_pufb2d_star.mat'),'U','alpha');
    assert(norm(methods(2).U-old.U,'fro')<1e-12);
    historical=load(fullfile(root,'pufb_otfs_v2','results','r2_2c_receiver_aware_fixed_budget','r22c_pufb_rx_star.mat'),'U');
    assert(norm(methods(3).U-historical.U,'fro')<1e-12);
    save(fullfile(output,'r22d_supports.mat'),'methods','-v7.3');
    write_arm_freeze(output,'r22d_swap_control_freeze.json',methods(6),[]);
    for axis={'nu','tau'}
        meta=axis_method.(upper_initial(axis{1}));
        T=table((0:3)',meta.h(:),meta.g(:),'VariableNames',{'Tap','Lowpass','Highpass'});
        writetable(T,fullfile(output,['r22d_filter_',axis{1},'.csv']));
    end
    write_status(output,'VAL_SELECTED_BOTH_ARMS_AND_TRAIN_COMPARATORS_FROZEN',false,false);
    result=struct('shared_alpha',shared_method.AlphaNu,'axis_alpha_nu',axis_method.AlphaNu,'axis_alpha_tau',axis_method.AlphaTau, ...
        'shared_val_jrx',shared_val_jrx,'axis_val_jrx',axis_val_jrx,'axis_global_candidates',height(global_surface),'axis_refined_candidates',height(refined));
end

function name=upper_initial(value)
    name=[upper(value(1)),value(2:end)];
end

function T=decorate(T,label)
    T.Stage=repmat(label,height(T),1);
end

function path=candidate_file(output,prefix,index)
    path=fullfile(output,prefix,sprintf('candidate_%04d.mat',index));
end

function T=evaluate_pairs(pairs,bank_values,singular,baseline,snr_db,p,output,prefix)
    directory=fullfile(output,prefix);if ~isfolder(directory),mkdir(directory);end
    rows=cell(size(pairs,1),7);
    for index=1:size(pairs,1)
        path=candidate_file(output,prefix,index);
        if isfile(path)
            loaded=load(path,'candidate');candidate=loaded.candidate;assert(max(abs(candidate.Pair-pairs(index,:)))<1e-13);
        else
            candidate=r22d_model('score',pairs(index,:),bank_values,singular,baseline,snr_db,p,[]);
            save(path,'candidate','-v7.3');
        end
        rows(index,:)={index,pairs(index,1),pairs(index,2),candidate.Jrx,candidate.Ratios(1),candidate.Ratios(2),candidate.Ratios(3)};
        fprintf('R22D %s %d/%d nu=%.9f tau=%.9f J=%.9g\n',prefix,index,size(pairs,1),pairs(index,1),pairs(index,2),candidate.Jrx);
    end
    T=cell2table(rows,'VariableNames',{'CandidateID','AlphaNu','AlphaTau','Jrx','Low','Mid','High'});
end

function centers=shared_neighborhoods(T,p)
    local=false(height(T),1);
    for index=1:height(T),neighbors=[mod(index-2,height(T))+1,mod(index,height(T))+1];local(index)=all(T.Jrx(index)<T.Jrx(neighbors));end
    ranked=[sortrows(T(local,:),{'Jrx','AlphaNu'});sortrows(T(~local,:),{'Jrx','AlphaNu'})];centers=[];
    for index=1:height(ranked)
        value=ranked.AlphaNu(index);distance=abs(value-centers);distance=min(distance,pi-distance);
        if isempty(centers)||all(distance>=3*pi/61-1e-12),centers(end+1)=value;end %#ok<AGROW>
        if numel(centers)==p.shared_neighborhoods,break;end
    end
    assert(numel(centers)==3);
end

function centers=axis_neighborhoods(T,p)
    local=false(height(T),1);spacing=pi/17;
    for index=1:height(T)
        distance=abs([T.AlphaNu-T.AlphaNu(index),T.AlphaTau-T.AlphaTau(index)]);distance=min(distance,pi-distance);
        neighbors=max(distance,[],2)<=spacing+1e-12&max(distance,[],2)>1e-12;
        local(index)=all(T.Jrx(index)<=T.Jrx(neighbors));
    end
    ranked=[sortrows(T(local,:),{'Jrx','AlphaNu','AlphaTau'});sortrows(T(~local,:),{'Jrx','AlphaNu','AlphaTau'})];centers=zeros(0,2);
    for index=1:height(ranked)
        pair=[ranked.AlphaNu(index),ranked.AlphaTau(index)];distance=abs(centers-pair);distance=min(distance,pi-distance);
        if isempty(centers)||all(max(distance,[],2)>=2*spacing-1e-12),centers(end+1,:)=pair;end %#ok<AGROW>
        if size(centers,1)==p.axis_neighborhoods,break;end
    end
    assert(size(centers,1)==4);
end

function values=bank(split,channels_per_velocity,p,physical,output)
    manifest_path=fullfile(output,['r22d_',split,'_manifest.csv']);T=r22d_model('read',manifest_path);
    selected=T.ChannelID<=channels_per_velocity;Tselected=T(selected,:);values=complex(zeros(288,288,height(Tselected)));
    bank_dir=fullfile(output,[split,'_channels']);if ~isfolder(bank_dir),mkdir(bank_dir);end
    for index=1:height(Tselected)
        velocity=Tselected.VelocityIndex(index);channel=Tselected.ChannelID(index);path=fullfile(bank_dir,sprintf('v%d_c%03d.mat',velocity,channel));
        if isfile(path),loaded=load(path,'G','seed');G=loaded.G;assert(loaded.seed==Tselected.Seed(index));
        else
            seed=Tselected.Seed(index);G=generate_channel_r20a(24,12,p.velocity_mps(velocity),seed,physical);save(path,'G','seed','-v7.3');
        end
        values(:,:,index)=G;T.ChannelGenerated(T.VelocityIndex==velocity&T.ChannelID==channel)=true;
    end
    writetable(T,manifest_path);
end

function singular=singular_bank(values)
    singular=zeros(288,size(values,3));
    for index=1:size(values,3),singular(:,index)=svd(values(:,:,index));end
end

function write_arm_freeze(output,name,method,val_jrx)
    r22a_io('json',fullfile(output,name),struct('Method',method.Name,'AlphaNu',method.AlphaNu,'AlphaTau',method.AlphaTau, ...
        'ValJrx',val_jrx,'K',nnz(method.Mask),'FilterLength',4,'Jnu',2,'Jtau',3,'StructuralResiduals',method.Residuals, ...
        'support_derived_from','R22D_FULL_TRAIN_ONLY','TEST_GENERATED',false));
end

function receipt=freeze_test(audit,output,p)
    path=fullfile(output,'r22d_temporal_freeze_receipt.json');assert(~isfile(path),'R22D:FreezeExists','Refusing to replace TEST freeze.');
    assert(isempty(dir(fullfile(output,'test_channels','*.mat'))));
    loaded=load(fullfile(output,'r22d_supports.mat'),'methods');methods=loaded.methods;
    for index=1:6
        checked=r22d_model('transform',methods(index).AlphaNu,methods(index).AlphaTau);
        assert(isequal(checked.U,methods(index).U)&&nnz(methods(index).Mask)==p.K);
    end
    files=dir(fullfile(audit,'*.m'));
    for index=1:numel(files),assert(isempty(checkcode(fullfile(audit,files(index).name),'-id')));end
    r22a_io('verify',fullfile(output,'r22d_code_hashes.csv'));
    r22a_io('verify',fullfile(output,'r22d_history_preservation_receipt.csv'));
    paths={fullfile(output,'r22d_supports.mat');fullfile(output,'r22d_primary_arms.mat'); ...
        fullfile(output,'r22d_shared_freeze.json');fullfile(output,'r22d_axis_specific_freeze.json'); ...
        fullfile(output,'r22d_swap_control_freeze.json');fullfile(output,'r22d_code_hashes.csv'); ...
        fullfile(output,'r22d_preregistration.json');fullfile(output,'r22d_train_candidates_frozen.mat'); ...
        fullfile(output,'r22d_val_selection.csv');fullfile(output,'r22d_history_preservation_receipt.csv')};
    hashes=r22d_model('hashes',paths);writetable(hashes,fullfile(output,'r22d_test_input_hashes.csv'));
    receipt=struct('CANDIDATES_FROZEN',true,'SUPPORTS_FROZEN',true,'COMPARATORS_FROZEN',true,'CODE_FROZEN',true, ...
        'STATISTICS_FROZEN',true,'FIGURES_REPORTS_FROZEN',true,'TEST_GENERATED',false,'TEST_ACCESSED',false, ...
        'utc',char(datetime('now','TimeZone','UTC')),'input_hashes_sha256',r22a_io('hash',fullfile(output,'r22d_test_input_hashes.csv')), ...
        'test_manifest_sha256',r22a_io('hash',fullfile(output,'r22d_test_manifest.csv')));
    r22a_io('json',path,receipt);persisted=jsondecode(fileread(path));assert(isequal(receipt,persisted));
    write_status(output,'COMPLETE_FREEZE_BEFORE_TEST',false,false);
end

function write_status(output,label,test_generated,test_accessed)
    r22a_io('json',fullfile(output,'r22d_execution_status.json'),struct('status',label, ...
        'VAL_SELECTION_ACCESS_COUNT',double(~any(strcmp(label,{'PREREGISTERED_BEFORE_NEW_BANKS'}))), ...
        'TEST_GENERATED',test_generated,'TEST_ACCESSED',test_accessed,'TEST_USED_FOR_SELECTION',false,'BER_RUN',false,'PCG_RUN',false));
end