function result=run_r22dcb(root,stage)
    pkg=fileparts(mfilename('fullpath'));pkg=fileparts(pkg);
    audit=fullfile(pkg,'04_ber_confirmation');
    output=fullfile(root,'pufb_otfs_v2','results','r2_2d_confirm_ber');
    source=fullfile(root,'pufb_otfs_v2','results','r2_2d_axis_specific');
    addpath(audit,fullfile(pkg,'02_shared_lib'), ...
        fullfile(root,'pufb_otfs_v2','audit','r2_2c_confirm_strict'),fullfile(pkg,'03_axis_specific_pufb'), ...
        fullfile(pkg,'01_channel_generator'));
    if ~isfolder(output),mkdir(output);end
    if strcmp(stage,'PREPARE')
        result=prepare(root,audit,output,source);return
    end
    p=jsondecode(fileread(fullfile(output,'r22dcb_preregistration.json')));
    receipt_path=fullfile(output,'r22dcb_temporal_freeze_receipt.json');receipt=jsondecode(fileread(receipt_path));
    assert(receipt.METHODS_FROZEN&&receipt.ALPHAS_FROZEN&&receipt.SUPPORTS_FROZEN&&receipt.CODE_FROZEN&&receipt.BER_RULES_FROZEN ...
        &&~receipt.PHASE_A_CHANNELS_GENERATED&&~receipt.PHASE_B_CHANNELS_GENERATED);
    assert(strcmp(receipt.code_hashes_sha256,r22a_io('hash',fullfile(output,'r22dcb_code_hashes.csv'))));
    assert(strcmp(receipt.input_hashes_sha256,r22a_io('hash',fullfile(output,'r22dcb_input_hashes.csv'))));
    verify(output);loaded=load(fullfile(output,'r22dcb_frozen_methods.mat'),'methods');methods=loaded.methods;
    physical=jsondecode(fileread(fullfile(root,'pufb_otfs_v2','audit','r2_0a_abid_true2d_pilot','r20a_preregistration.json')));
    if strcmp(stage,'A')
        result=phase_a(methods,p,physical,output);verify(output);
        manifest=r22a_io('manifest',output,audit);fprintf('PHASE A %s; %d artifacts\n',result.R22D_CONFIRMATION,height(manifest));return
    end
    assert(strcmp(stage,'B'),'R22DCB:Stage','Use PREPARE, A or B.');
    decision=jsondecode(fileread(fullfile(output,'r22dcb_phase_a_decision.json')));
    assert(strcmp(decision.R22D_CONFIRMATION,'CONFIRMED'),'R22DCB:Authorization','PHASE A did not confirm; BER forbidden.');
    authorization_path=fullfile(output,'r22dcb_phase_b_authorization.json');
    if ~isfile(authorization_path)
        assert(isempty(dir(fullfile(output,'phase_b_channels','*.mat'))));verify(output);
        authorization=struct('R22D_CONFIRMATION','CONFIRMED','METHODS_UNCHANGED_AFTER_PHASE_A',true, ...
            'SUPPORTS_UNCHANGED_AFTER_PHASE_A',true,'BER_RULES_UNCHANGED_AFTER_PHASE_A',true, ...
            'PHASE_B_CHANNELS_GENERATED',false,'utc',char(datetime('now','TimeZone','UTC')), ...
            'phase_a_decision_sha256',r22a_io('hash',fullfile(output,'r22dcb_phase_a_decision.json')), ...
            'phase_a_manifest_sha256',r22a_io('hash',fullfile(output,'r22dcb_phase_a_manifest.csv')));
        r22a_io('json',authorization_path,authorization);
    end
    authorization=jsondecode(fileread(authorization_path));
    assert(strcmp(authorization.phase_a_decision_sha256,r22a_io('hash',fullfile(output,'r22dcb_phase_a_decision.json'))) ...
        &&authorization.METHODS_UNCHANGED_AFTER_PHASE_A&&authorization.SUPPORTS_UNCHANGED_AFTER_PHASE_A ...
        &&authorization.BER_RULES_UNCHANGED_AFTER_PHASE_A&&~authorization.PHASE_B_CHANNELS_GENERATED);
    verify(output);write_status(output,'PHASE_B_AUTHORIZED_AND_RUNNING',true,true);
    r22dcb_ber(methods,p,physical,output);result=r22dcb_finish(methods,p,output);
    verify(output);write_status(output,'COMPLETED',true,true);
    manifest=r22a_io('manifest',output,audit);r22a_io('verify',fullfile(output,'artifact_manifest.csv'));
    fprintf('R22DCB COMPLETE %s; %d artifacts\n',result.PRIMARY_INTERPRETATION,height(manifest));
end

function p=prepare(root,audit,output,source)
    assert(isempty(dir(fullfile(output,'phase_a_channels','*.mat')))&&isempty(dir(fullfile(output,'phase_b_channels','*.mat'))));
    assert(~isfile(fullfile(output,'r22dcb_temporal_freeze_receipt.json')),'R22DCB:AlreadyFrozen','Refusing to overwrite complete freeze.');
    r22a_io('verify',fullfile(source,'artifact_manifest.csv'));
    r22a_io('verify',fullfile(source,'r22d_history_preservation_receipt.csv'));
    history=r22dcb_math('read',fullfile(source,'artifact_manifest.csv'));history=history(:,{'AbsolutePath','SHA256'});
    earlier=r22dcb_math('read',fullfile(source,'r22d_history_preservation_receipt.csv'));
    history=[history;earlier(:,{'AbsolutePath','SHA256'})];writetable(history,fullfile(output,'r22dcb_history_hashes.csv'));
    loaded=load(fullfile(source,'r22d_supports.mat'),'methods');methods=loaded.methods([1,3,5,2]);
    names=["D4_2D_PARTIAL_MALLAT","R2.2C_SHARED_RX","R2.2D_AXIS_RX","R2.2A_PUFB2D_STAR"];
    assert(methods(3).AlphaNu==0&&abs(methods(3).AlphaTau-.78539816339745)<1e-14 ...
        &&abs(methods(2).AlphaNu-.78797323934301)<1e-14&&methods(2).AlphaNu==methods(2).AlphaTau);
    method_metadata=cell(4,1);support_metadata=cell(4,1);
    for index=1:4
        methods(index).Name=names(index);original=methods(index);
        verification=r22d_model('transform',original.AlphaNu,original.AlphaTau);
        assert(norm(verification.U-original.U,'fro')<1e-12&&max(original.Residuals)<1e-11&&nnz(original.Mask)==33178);
        assert(isequal(verification.Pnu,original.Pnu)&&isequal(verification.Ptau,original.Ptau) ...
            &&isequal(verification.Nu.h,original.Nu.h)&&isequal(verification.Tau.h,original.Tau.h));
        method_metadata{index}=struct('Method',names(index),'AlphaNu',original.AlphaNu,'AlphaTau',original.AlphaTau, ...
            'U_SHA256',r22dcb_math('arrayhash',original.U),'Pnu_SHA256',r22dcb_math('arrayhash',original.Pnu), ...
            'Ptau_SHA256',r22dcb_math('arrayhash',original.Ptau), ...
            'FilterNu_SHA256',r22dcb_math('arrayhash',[original.Nu.h,original.Nu.g]), ...
            'FilterTau_SHA256',r22dcb_math('arrayhash',[original.Tau.h,original.Tau.g]), ...
            'PolyphaseNu_SHA256',r22dcb_math('arrayhash',original.Nu.E), ...
            'PolyphaseTau_SHA256',r22dcb_math('arrayhash',original.Tau.E),'Residuals',original.Residuals);
        support_metadata{index}=struct('Method',names(index),'K',33178, ...
            'Mask_SHA256',r22dcb_math('arrayhash',uint8(original.Mask(:))), ...
            'SourceMAT_SHA256',r22a_io('hash',fullfile(source,'r22d_supports.mat')));
    end
    save(fullfile(output,'r22dcb_frozen_methods.mat'),'methods','-v7.3');
    r22a_io('json',fullfile(output,'r22dcb_method_freeze.json'),vertcat(method_metadata{:}));
    r22a_io('json',fullfile(output,'r22dcb_support_freeze.json'),vertcat(support_metadata{:}));
    p=struct('campaign','R2.2D-CONFIRM_PLUS_BER','physical_scenario','ABID_EVA','M',24,'N',12,'L',288, ...
        'K',33178,'Lw',4,'Jnu',2,'Jtau',3,'phase_a_per_velocity',200,'phase_b_per_velocity',100, ...
        'phase_a_snr_db',3:3:24,'ber_snr_db',0:1:24, ...
        'velocity_mps',[50/3.6,100/3.6,500/3.6,340.3,1701.5], ...
        'velocity_labels',{{'50 km/h','100 km/h','500 km/h','Mach 1','Mach 5'}}, ...
        'seeds',struct('A',132100000,'B',132200000,'DATA',133100000,'NOISE',134100000,'MSE_BOOT',135100000,'BER_BOOT',135200000,'FULL_DATA',136100000,'FULL_NOISE',137100000), ...
        'physical_seed_rule','base+1000*velocity+local_channel; A and B completely disjoint', ...
        'frame_seed_rule','base+100000*velocity+1000*SNRdB+local_channel; RandStream mt19937ar; reset per channel/condition', ...
        'modulation','Gray QPSK (1-2*real_bit+j*(1-2*imag_bit))/sqrt(2); bits/sign mapping unchanged', ...
        'bootstrap_replicates',10000,'phase_a_weighting','Equal LOW/MID/HIGH channel averages then equal5velocity averages', ...
        'phase_a_bootstrap','Physical paired channels; independently resample200 within each velocity; never SNR/frames independent units', ...
        'phase_a_success','Primary axis-minus-shared upper95percentile<0. Mean>=0 NOT_CONFIRMED, otherwise unresolvedCI INCONCLUSIVE. No secondary override.', ...
        'ber_target_errors',300,'ber_cap_bits',50000000,'ber_batch_rounds',16,'ber_frame_bits',576, ...
        'ber_stop_rule','COMMON rounds for all4methods; stop iff ALL3primary aggregate errors>=300 OR floor(cap/(576*100)) rounds reached; secondary never controls stop', ...
        'ber_aggregate_weighting','Average per-velocity ratios errors/bits with equal velocity weights; high-mobility equal Mach1/Mach5', ...
        'ber_targets',[1e-2,1e-3,1e-4],'primary_target',1e-3,'secondary_target',1e-4, ...
        'crossing_rule','First descending adjacent pair of measured0:24dB points bracketing target with positive BER and each >=300 measured aggregate errors; interpolate SNR versus log10BER only. No smoothing, spline, extrapolation or floor replacement.', ...
        'bootstrap_crossing_rule','Recompute count ratios per resampled physical channel, per-velocity then group mean curves; same resamples across methods/SNR; original measured error>=300 eligibility fixed; positive bootstrap points required', ...
        'bootstrap_resolution_rule','Require original paired crossings AND at least95percent of10000 resamples resolved; otherwise NOT_RESOLVED. PercentileCI conditional on resolved samples, resolution count/fraction reported explicitly.', ...
        'gain_success','axis-vs-shared1e-3 lower95percentile>0; no0.5/1/2dB success thresholds', ...
        'descriptive_gain_bins','<0.25 VERY_SMALL;[0.25,0.5) SMALL;[0.5,1) MODERATE;[1,2] LARGE;>2 VERY_LARGE; negative NONPOSITIVE', ...
        'high_mobility_group',[4,5],'high_mobility_dominance_rule','Aggregate1e-3 supported AND highmobility supported AND pooled lowmobility(1:3) not supported with resolved mean lower than highmobility gain', ...
        'curve_ci_rule','Per-velocity binomialWilson; group interval average of Bonferroni simultaneous per-velocity Wilson bounds (family alpha.05); unconditional method comparison uses channel bootstrap.', ...
        'full_control_channels_per_velocity',10,'full_control_snr_db',[6,15,24],'full_control_rounds',64, ...
        'full_control_pass_rule','Aggregate equal velocity and SNR paired full-axis minus full-shared BER bootstrapCI includes zero. Diagnostic only; unitary trace-MSE invariance does not guarantee finite-QPSK BER invariance.', ...
        'context_snr_db',3:3:24,'full_mse_tolerance',1e-12,'NO_OPTIMIZATION',true,'MANUSCRIPT_MODIFIED',false);
    prereg_path=fullfile(output,'r22dcb_preregistration.json');
    if isfile(prereg_path),assert(isequaln(jsondecode(fileread(prereg_path)),jsondecode(jsonencode(p))));
    else,r22a_io('json',prereg_path,p);end
    phase=repelem(["A";"B"],[1000;500]);velocity=[repelem((1:5)',200);repelem((1:5)',100)];
    channel=[repmat((1:200)',5,1);repmat((1:100)',5,1)];seed=[p.seeds.A+1000*velocity(1:1000)+channel(1:1000); ...
        p.seeds.B+1000*velocity(1001:end)+channel(1001:end)];
    existing=dir(fullfile(root,'pufb_otfs_v2','results','**','*manifest.csv'));
    for index=1:numel(existing)
        path=fullfile(existing(index).folder,existing(index).name);if startsWith(path,[output,filesep]),continue;end
        T=r22dcb_math('read',path);if ismember('Seed',T.Properties.VariableNames),assert(isempty(intersect(seed,T.Seed)));end
    end
    manifest=table(phase,velocity,channel,seed,false(1500,1),'VariableNames',{'Phase','VelocityIndex','ChannelID','Seed','ChannelGenerated'});
    writetable(manifest,fullfile(output,'r22dcb_seed_manifest.csv'));
    writetable(manifest(1:1000,:),fullfile(output,'r22dcb_phase_a_manifest.csv'));
    r22a_io('json',fullfile(output,'r22dcb_preflight_math.json'),r22dcb_math('selftest'));
    r22dcb_finish('SELFTEST',p,output);
    files=dir(fullfile(audit,'*.m'));paths=cell(numel(files),1);
    for index=1:numel(files)
        paths{index}=fullfile(audit,files(index).name);d=checkcode(paths{index},'-id');
        for diagnostic=1:numel(d),fprintf('%s:%d %s\n',files(index).name,d(diagnostic).line,d(diagnostic).message);end
        assert(isempty(d),'R22DCB:Checkcode','Pre-generation source diagnostic.');
    end
    paths=[paths;{which('r22a_io');which('fir_mallat_r22a');which('r22d_model');which('confirm_exact_mse');which('generate_channel_r20a'); ...
        fullfile(root,'pufb_otfs_v2','audit','r2_0a_abid_true2d_pilot','r20a_preregistration.json')}];
    hashes=r22d_model('hashes',paths);writetable(hashes,fullfile(output,'r22dcb_code_hashes.csv'));
    inputs={fullfile(output,'r22dcb_frozen_methods.mat');fullfile(output,'r22dcb_method_freeze.json');fullfile(output,'r22dcb_support_freeze.json'); ...
        fullfile(output,'r22dcb_preregistration.json');fullfile(output,'r22dcb_seed_manifest.csv');fullfile(output,'r22dcb_history_hashes.csv')};
    writetable(r22d_model('hashes',inputs),fullfile(output,'r22dcb_input_hashes.csv'));
    receipt=struct('METHODS_FROZEN',true,'ALPHAS_FROZEN',true,'SUPPORTS_FROZEN',true,'CODE_FROZEN',true,'BER_RULES_FROZEN',true, ...
        'PHASE_A_CHANNELS_GENERATED',false,'PHASE_A_CHANNELS_ACCESSED',false,'PHASE_B_CHANNELS_GENERATED',false,'PHASE_B_CHANNELS_ACCESSED',false, ...
        'utc',char(datetime('now','TimeZone','UTC')),'code_hashes_sha256',r22a_io('hash',fullfile(output,'r22dcb_code_hashes.csv')), ...
        'input_hashes_sha256',r22a_io('hash',fullfile(output,'r22dcb_input_hashes.csv')));
    receipt_path=fullfile(output,'r22dcb_temporal_freeze_receipt.json');r22a_io('json',receipt_path,receipt);
    assert(isequal(receipt,jsondecode(fileread(receipt_path))));verify(output);write_status(output,'ALL_PHASE_A_AND_B_RULES_FROZEN_BEFORE_CHANNELS',false,false);
end

function verify(output)
    for name={'r22dcb_history_hashes.csv','r22dcb_code_hashes.csv','r22dcb_input_hashes.csv'},r22a_io('verify',fullfile(output,name{1}));end
end

function decision=phase_a(methods,p,physical,output)
    directory=fullfile(output,'phase_a_channels');if ~isfolder(directory),mkdir(directory);end
    manifest=r22dcb_math('read',fullfile(output,'r22dcb_phase_a_manifest.csv'));
    assert(isequal(manifest.Seed,p.seeds.A+1000*repelem((1:5)',200)+repmat((1:200)',5,1)));
    rows=cell(32000,8);full_rows=cell(32000,5);frob_rows=cell(4000,5);metric_cursor=0;channel_cursor=0;
    write_status(output,'PHASE_A_RUNNING',true,false);
    for index=1:1000
        path=fullfile(directory,sprintf('channel_%04d.mat',index));
        if isfile(path),loaded=load(path,'record','seed');record=loaded.record;assert(loaded.seed==manifest.Seed(index));
        else
            seed=manifest.Seed(index);velocity=manifest.VelocityIndex(index);G=generate_channel_r20a(24,12,p.velocity_mps(velocity),seed,physical);
            record=struct('SparseMSE',zeros(4,8),'FullMSE',zeros(4,8),'Frobenius',zeros(4,1));
            for method=1:4
                H=methods(method).U'*G*methods(method).U;assert(nnz(methods(method).Mask)==33178);
                [record.SparseMSE(method,:),record.FullMSE(method,:)]=confirm_exact_mse(H,methods(method).Mask,p.phase_a_snr_db);
                record.Frobenius(method)=norm(H-H.*methods(method).Mask,'fro')/norm(H,'fro');
            end
            assert(max(abs(record.FullMSE-record.FullMSE(1,:)),[],'all')<=p.full_mse_tolerance);save(path,'G','record','seed','-v7.3');
        end
        for method=1:4
            name=methods(method).Name;velocity=manifest.VelocityIndex(index);channel=manifest.ChannelID(index);
            channel_cursor=channel_cursor+1;frob_rows(channel_cursor,:)={name,velocity,channel,p.K,record.Frobenius(method)};
            for point=1:8
                metric_cursor=metric_cursor+1;rows(metric_cursor,:)={name,velocity,channel,p.phase_a_snr_db(point),p.K, ...
                    record.SparseMSE(method,point),record.FullMSE(method,point),record.SparseMSE(method,point)-record.FullMSE(method,point)};
                full_rows(metric_cursor,:)={name,velocity,channel,p.phase_a_snr_db(point),record.FullMSE(method,point)};
            end
        end
        manifest.ChannelGenerated(index)=true;
        if mod(index,10)==0||index==1000,writetable(manifest,fullfile(output,'r22dcb_phase_a_manifest.csv'));fprintf('R22DCB phaseA %d/1000\n',index);end
    end
    mse=cell2table(rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','RetainedEntries','MSESparse','MSEFull','DeltaMSE'});
    full=cell2table(full_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','FullMSE'});
    frob=cell2table(frob_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','RetainedEntries','FrobeniusError'});
    writetable(mse,fullfile(output,'r22dcb_phase_a_exact_mse.csv'));writetable(full,fullfile(output,'r22dcb_phase_a_full_mse.csv'));
    writetable(frob,fullfile(output,'r22dcb_phase_a_frobenius.csv'));
    decision=r22dcb_finish('A',p,output);write_status(output,['PHASE_A_',decision.R22D_CONFIRMATION],true,false);
end

function write_status(output,label,phase_a,phase_b)
    r22a_io('json',fullfile(output,'r22dcb_execution_status.json'),struct('status',label, ...
        'PHASE_A_CHANNELS_GENERATED',phase_a,'PHASE_A_CHANNELS_ACCESSED',phase_a, ...
        'PHASE_B_CHANNELS_GENERATED',phase_b,'PHASE_B_CHANNELS_ACCESSED',phase_b,'TEST_USED_FOR_SELECTION',false, ...
        'NEW_TRANSFORM_TRAINED',false,'SUPPORT_REBUILT',false,'MANUSCRIPT_MODIFIED',false));
end