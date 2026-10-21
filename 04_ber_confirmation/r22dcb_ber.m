function r22dcb_ber(methods,p,physical,output)
    directory=fullfile(output,'phase_b_channels');if ~isfolder(directory),mkdir(directory);end
    seed_manifest=r22dcb_math('read',fullfile(output,'r22dcb_seed_manifest.csv'));
    manifest=seed_manifest(seed_manifest.Phase=="B",:);assert(height(manifest)==500);
    manifest_path=fullfile(output,'r22dcb_ber_manifest.csv');
    if isfile(manifest_path),manifest=r22dcb_math('read',manifest_path);end
    assert(isequal(manifest.Seed,p.seeds.B+1000*repelem((1:5)',100)+repmat((1:100)',5,1)));
    counts=table();per_channel=table();context=table();
    for velocity=1:5
        physical_channels=cell(100,1);effective=cell(4,100);context_rows=cell(400,8);context_cursor=0;
        context_path=fullfile(output,sprintf('ber_context_v%d.mat',velocity));
        if isfile(context_path),loaded=load(context_path,'context_table');context_table=loaded.context_table;else,context_table=[];end
        for channel=1:100
            index=(velocity-1)*100+channel;path=fullfile(directory,sprintf('channel_%04d.mat',index));
            if isfile(path),loaded=load(path,'G','seed');G=loaded.G;assert(loaded.seed==manifest.Seed(index));
            else,seed=manifest.Seed(index);G=generate_channel_r20a(24,12,p.velocity_mps(velocity),seed,physical);save(path,'G','seed','-v7.3');end
            physical_channels{channel}=G;manifest.ChannelGenerated(index)=true;
            for method=1:4
                U=methods(method).U;H=U'*G*U;effective{method,channel}=H;
                if isempty(context_table)
                    [sparse_mse,full_mse]=confirm_exact_mse(H,methods(method).Mask,p.context_snr_db);
                    delta=sparse_mse-full_mse;balanced=(mean(delta(1:3))+mean(delta(4:6))+mean(delta(7:8)))/3;
                    energy=abs(H(:)).^2;probability=energy/sum(energy);cumulative=cumsum(sort(probability,'descend'));
                    context_cursor=context_cursor+1;context_rows(context_cursor,:)={methods(method).Name,velocity,channel,p.K,balanced, ...
                        norm(H-H.*methods(method).Mask,'fro')/norm(H,'fro'),find(cumulative>=.99,1),sum(probability.^2)};
                end
            end
        end
        if isempty(context_table)
            context_table=cell2table(context_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','RetainedEntries','MeanDeltaMSE','FrobeniusError','K99','J4'});
            save(context_path,'context_table');
        end
        context=[context;context_table]; %#ok<AGROW>
        writetable(context,fullfile(output,'r22dcb_ber_context_metrics.csv'));writetable(manifest,manifest_path);
        for snr=p.ber_snr_db(:).'
            checkpoint=fullfile(output,sprintf('ber_v%d_snr%02d.mat',velocity,snr));
            if isfile(checkpoint)
                loaded=load(checkpoint,'summary_table','channel_table');summary_table=loaded.summary_table;channel_table=loaded.channel_table;
            else
                variance=10^(-snr/10);signal_map=cell(4,100);noise_map=signal_map;
                for channel=1:100
                    for method=1:4
                        H=effective{method,channel};Hs=H.*methods(method).Mask;receiver=r22dcb_math('receiver',Hs,variance);
                        signal_map{method,channel}=receiver*H;noise_map{method,channel}=receiver*methods(method).U';
                    end
                end
                data=cell(100,1);noise=data;
                for channel=1:100
                    data{channel}=RandStream('mt19937ar','Seed',p.seeds.DATA+100000*velocity+1000*snr+channel);
                    noise{channel}=RandStream('mt19937ar','Seed',p.seeds.NOISE+100000*velocity+1000*snr+channel);
                end
                errors=zeros(4,100);rounds=0;maximum=floor(p.ber_cap_bits/(p.ber_frame_bits*100));
                while rounds<maximum&&~all(sum(errors(1:3,:),2)>=p.ber_target_errors)
                    batch=min(p.ber_batch_rounds,maximum-rounds);rounds=rounds+batch;
                    for channel=1:100
                        bits=logical(randi(data{channel},[0,1],576,batch));real_bits=bits(1:2:end,:);imag_bits=bits(2:2:end,:);
                        symbols=(1-2*double(real_bits)+1i*(1-2*double(imag_bits)))/sqrt(2);
                        physical_noise=sqrt(variance/2)*(randn(noise{channel},288,batch)+1i*randn(noise{channel},288,batch));
                        for method=1:4
                            estimate=signal_map{method,channel}*symbols+noise_map{method,channel}*physical_noise;
                            errors(method,channel)=errors(method,channel)+nnz((real(estimate)<0)~=real_bits)+nnz((imag(estimate)<0)~=imag_bits);
                        end
                    end
                end
                summary_rows=cell(4,13);channel_rows=cell(400,9);cursor=0;
                for method=1:4
                    total_errors=sum(errors(method,:));bits=p.ber_frame_bits*rounds*100;interval=r22dcb_math('wilson',total_errors,bits);
                    summary_rows(method,:)={methods(method).Name,velocity,snr,total_errors,bits,rounds,100,total_errors==0, ...
                        rounds==maximum,total_errors/bits,interval(1),interval(2),p.K};
                    for channel=1:100
                        cursor=cursor+1;channel_rows(cursor,:)={methods(method).Name,velocity,channel,snr,errors(method,channel), ...
                            p.ber_frame_bits*rounds,rounds,p.K,manifest.Seed((velocity-1)*100+channel)};
                    end
                end
                summary_table=cell2table(summary_rows,'VariableNames',{'Method','VelocityIndex','SNRdB','Errors','Bits','FramesPerChannel', ...
                    'ChannelsUsed','ZeroError','CapReached','BER','WilsonLower95','WilsonUpper95','RetainedEntries'});
                channel_table=cell2table(channel_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','Errors','Bits','Frames','RetainedEntries','Seed'});
                save(checkpoint,'summary_table','channel_table','-v7.3');
            end
            counts=[counts;summary_table];per_channel=[per_channel;channel_table]; %#ok<AGROW>
            writetable(counts,fullfile(output,'r22dcb_ber_counts.csv'));writetable(per_channel,fullfile(output,'r22dcb_ber_per_channel.csv'));
            fprintf('R22DCB BER velocity=%d SNR=%d rounds=%d errors=%s\n',velocity,snr,summary_table.FramesPerChannel(1),mat2str(summary_table.Errors.'));
        end
        full_control(methods,p,physical_channels,velocity,output);
    end
    clear signal_map noise_map effective physical_channels
end

function full_control(methods,p,physical_channels,velocity,output)
    path=fullfile(output,sprintf('full_control_v%d.mat',velocity));
    if isfile(path),return;end
    rows=cell(2*p.full_control_channels_per_velocity*numel(p.full_control_snr_db),8);cursor=0;
    for snr=p.full_control_snr_db(:).'
        variance=10^(-snr/10);
        for channel=1:p.full_control_channels_per_velocity
            G=physical_channels{channel};data=RandStream('mt19937ar','Seed',p.seeds.FULL_DATA+100000*velocity+1000*snr+channel);
            noise=RandStream('mt19937ar','Seed',p.seeds.FULL_NOISE+100000*velocity+1000*snr+channel);
            bits=logical(randi(data,[0,1],576,p.full_control_rounds));real_bits=bits(1:2:end,:);imag_bits=bits(2:2:end,:);
            symbols=(1-2*double(real_bits)+1i*(1-2*double(imag_bits)))/sqrt(2);
            physical_noise=sqrt(variance/2)*(randn(noise,288,p.full_control_rounds)+1i*randn(noise,288,p.full_control_rounds));
            for method=[2,3]
                U=methods(method).U;H=U'*G*U;receiver=r22dcb_math('receiver',H,variance);
                estimate=receiver*(H*symbols+U'*physical_noise);errors=nnz((real(estimate)<0)~=real_bits)+nnz((imag(estimate)<0)~=imag_bits);
                cursor=cursor+1;rows(cursor,:)={methods(method).Name,velocity,channel,snr,errors,576*p.full_control_rounds, ...
                    p.full_control_rounds,errors/(576*p.full_control_rounds)};
            end
        end
    end
    full_table=cell2table(rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','Errors','Bits','Frames','BER'});
    save(path,'full_table');
end