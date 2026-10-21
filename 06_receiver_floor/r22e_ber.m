function r22e_ber(methods,etas,p,~,output)
    assert(isfile(fullfile(output,'r22e_pretest_freeze_receipt.json'))&&isfile(fullfile(output,'r22e_numerical_controls.json')));
    indices=[3,2,1,4,3];parameters=[etas(3),etas(2),etas(1),etas(4),0];names=["AXIS_REG","SHARED_REG","D4_REG","R22A_REG","AXIS_ETA0"];
    counts=table();per_channel=table();
    for velocity=1:5
        effective=cell(4,200);
        for channel=1:200
            loaded=load(fullfile(output,'test_channels',sprintf('channel_%04d.mat',(velocity-1)*200+channel)),'G');G=loaded.G;
            for method=1:4,effective{method,channel}=methods(method).U'*G*methods(method).U;end
        end
        for snr=p.ber_snr(:).'
            checkpoint=fullfile(output,sprintf('ber_v%d_snr%02d.mat',velocity,snr));
            if isfile(checkpoint),loaded=load(checkpoint,'summary_table','channel_table');summary_table=loaded.summary_table;channel_table=loaded.channel_table;
            else
                variance=10^(-snr/10);signal_map=cell(5,200);noise_map=signal_map;
                for channel=1:200
                    for variant=1:5
                        method=indices(variant);H=effective{method,channel};Hs=H.*methods(method).Mask;
                        W=r22e_math('receiver',Hs,variance,parameters(variant));signal_map{variant,channel}=W*H;noise_map{variant,channel}=W*methods(method).U';
                    end
                end
                data=cell(200,1);noise=data;
                for channel=1:200
                    data{channel}=RandStream('mt19937ar','Seed',p.seeds.DATA+100000*velocity+1000*snr+channel);
                    noise{channel}=RandStream('mt19937ar','Seed',p.seeds.NOISE+100000*velocity+1000*snr+channel);
                end
                errors=zeros(5,200);rounds=0;maximum=floor(p.ber_cap_bits/(576*200));
                while rounds<maximum&&~all(sum(errors([1,2,3,5],:),2)>=p.ber_target_errors)
                    batch=min(p.ber_batch_rounds,maximum-rounds);rounds=rounds+batch;
                    for channel=1:200
                        bits=logical(randi(data{channel},[0,1],576,batch));real_bits=bits(1:2:end,:);imag_bits=bits(2:2:end,:);
                        symbols=(1-2*double(real_bits)+1i*(1-2*double(imag_bits)))/sqrt(2);
                        physical_noise=sqrt(variance/2)*(randn(noise{channel},288,batch)+1i*randn(noise{channel},288,batch));
                        for variant=1:5
                            estimate=signal_map{variant,channel}*symbols+noise_map{variant,channel}*physical_noise;
                            errors(variant,channel)=errors(variant,channel)+nnz((real(estimate)<0)~=real_bits)+nnz((imag(estimate)<0)~=imag_bits);
                        end
                    end
                end
                rows=cell(5,14);channel_rows=cell(1000,9);cursor=0;
                for variant=1:5
                    error_count=sum(errors(variant,:));bit_count=576*200*rounds;interval=r22dcb_math('wilson',error_count,bit_count);
                    rows(variant,:)={names(variant),velocity,snr,error_count,bit_count,rounds,200,error_count==0,rounds==maximum, ...
                        error_count/bit_count,interval(1),interval(2),33178,parameters(variant)};
                    for channel=1:200
                        cursor=cursor+1;channel_rows(cursor,:)={names(variant),velocity,channel,snr,errors(variant,channel),576*rounds,rounds,33178,parameters(variant)};
                    end
                end
                summary_table=cell2table(rows,'VariableNames',{'Variant','VelocityIndex','SNRdB','Errors','Bits','FramesPerChannel','ChannelsUsed', ...
                    'ZeroError','CapReached','BER','Lower95','Upper95','RetainedEntries','Eta'});
                channel_table=cell2table(channel_rows,'VariableNames',{'Variant','VelocityIndex','ChannelID','SNRdB','Errors','Bits','Frames','RetainedEntries','Eta'});
                save(checkpoint,'summary_table','channel_table','-v7.3');
            end
            counts=[counts;summary_table];per_channel=[per_channel;channel_table]; %#ok<AGROW>
            writetable(counts,fullfile(output,'r22e_ber_counts.csv'));writetable(per_channel,fullfile(output,'r22e_ber_per_channel.csv'));
            fprintf('R22E BER velocity%d SNR%d rounds%d errors%s\n',velocity,snr,summary_table.FramesPerChannel(1),mat2str(summary_table.Errors.'));
        end
    end
end