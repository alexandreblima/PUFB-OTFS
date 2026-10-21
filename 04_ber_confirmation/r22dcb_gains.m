function result=r22dcb_gains(methods,p,per_channel,output)
    errors=zeros(100,25,4,5);bits=errors;
    for method=1:4
        for velocity=1:5
            selected=sortrows(per_channel(string(per_channel.Method)==methods(method).Name&per_channel.VelocityIndex==velocity,:),{'ChannelID','SNRdB'});
            assert(height(selected)==2500&&isequal(selected.ChannelID,repelem((1:100)',25)) ...
                &&isequal(selected.SNRdB,repmat((0:24)',100,1))&&all(selected.RetainedEntries==p.K));
            errors(:,:,method,velocity)=reshape(selected.Errors,25,100).';bits(:,:,method,velocity)=reshape(selected.Bits,25,100).';
        end
    end
    groups=struct('Name',{'AGGREGATE','HIGH_MOBILITY','LOW_MOBILITY'},'Velocities',{1:5,[4,5],1:3});
    for velocity=1:5,groups(end+1)=struct('Name',sprintf('VELOCITY_%d',velocity),'Velocities',velocity);end %#ok<AGROW>
    summary_rows=cell(numel(groups)*4*25,8);cross_rows=cell(numel(groups)*4*3,7);gain_rows=cell(numel(groups)*3*3,15);
    summary_cursor=0;cross_cursor=0;gain_cursor=0;
    for group=1:numel(groups)
        velocities=groups(group).Velocities;rate=zeros(4,25);eligible=false(4,25);lower=rate;upper=rate;
        boot=zeros(p.bootstrap_replicates,25,4);stream=RandStream('mt19937ar','Seed',p.seeds.BER_BOOT);
        for velocity=velocities
            draws=randi(stream,100,[100,p.bootstrap_replicates]);
            weights=accumarray([draws(:),repelem((1:p.bootstrap_replicates)',100)],1,[100,p.bootstrap_replicates]);
            for method=1:4
                current_errors=errors(:,:,method,velocity);current_bits=bits(:,:,method,velocity);
                rate(method,:)=rate(method,:)+sum(current_errors,1)./sum(current_bits,1)/numel(velocities);
                boot(:,:,method)=boot(:,:,method)+(weights'*current_errors)./(weights'*current_bits)/numel(velocities);
                for point=1:25
                    interval=r22dcb_math('wilson',sum(current_errors(:,point)),sum(current_bits(:,point)),.05/numel(velocities));
                    lower(method,point)=lower(method,point)+interval(1)/numel(velocities);
                    upper(method,point)=upper(method,point)+interval(2)/numel(velocities);
                end
            end
        end
        for method=1:4
            totals=sum(errors(:,:,method,velocities),[1,4]);totals=reshape(totals,1,25);eligible(method,:)=totals>=p.ber_target_errors;
            for point=1:25
                summary_cursor=summary_cursor+1;summary_rows(summary_cursor,:)={string(groups(group).Name),methods(method).Name,point-1, ...
                    rate(method,point),lower(method,point),upper(method,point),totals(point),eligible(method,point)};
            end
        end
        crossings=nan(4,3);bootstrap_crossings=nan(p.bootstrap_replicates,4,3);
        for method=1:4
            for target_index=1:3
                target=p.ber_targets(target_index);cross=r22dcb_math('crossing',0:24,rate(method,:),target,eligible(method,:));
                crossings(method,target_index)=cross.SNR;cross_cursor=cross_cursor+1;
                cross_rows(cross_cursor,:)={string(groups(group).Name),methods(method).Name,target,cross.SNR,string(cross.Status),cross.LeftIndex-1,cross.RightIndex-1};
                if isfinite(cross.SNR)
                    for replicate=1:p.bootstrap_replicates
                        sampled=r22dcb_math('crossing',0:24,boot(replicate,:,method),target,eligible(method,:));
                        bootstrap_crossings(replicate,method,target_index)=sampled.SNR;
                    end
                end
            end
        end
        for reference=[2,1,4]
            for target_index=1:3
                target=p.ber_targets(target_index);samples=bootstrap_crossings(:,reference,target_index)-bootstrap_crossings(:,3,target_index);
                valid=isfinite(samples);fraction=mean(valid);original=crossings(reference,target_index)-crossings(3,target_index);
                lower_gain=NaN;upper_gain=NaN;mean_gain=NaN;median_gain=NaN;status='NOT_RESOLVED';
                if isfinite(original)&&fraction>=.95
                    ordered=sort(samples(valid));mean_gain=mean(ordered);median_gain=median(ordered);
                    lower_gain=ordered(ceil(.025*numel(ordered)));upper_gain=ordered(ceil(.975*numel(ordered)));
                    status='NOT_SUPPORTED';if lower_gain>0,status='SUPPORTED';end
                end
                descriptive='NOT_RESOLVED';
                if isfinite(original)
                    if original<=0
                        descriptive='NONPOSITIVE';
                    elseif original<.25
                        descriptive='VERY_SMALL';
                    elseif original<.5
                        descriptive='SMALL';
                    elseif original<1
                        descriptive='MODERATE';
                    elseif original<=2
                        descriptive='LARGE';
                    else
                        descriptive='VERY_LARGE';
                    end
                end
                gain_cursor=gain_cursor+1;gain_rows(gain_cursor,:)={string(groups(group).Name),"AXIS_VS_"+methods(reference).Name,target, ...
                    crossings(reference,target_index),crossings(3,target_index),original,mean_gain,median_gain,lower_gain,upper_gain, ...
                    nnz(valid),p.bootstrap_replicates,fraction,string(status),string(descriptive)};
            end
        end
        fprintf('BER target crossing/bootstrap group %s done\n',groups(group).Name);
    end
    summary=cell2table(summary_rows,'VariableNames',{'Group','Method','SNRdB','BER','Lower95','Upper95','Errors','CrossingEligible'});
    crossings=cell2table(cross_rows,'VariableNames',{'Group','Method','TargetBER','CrossingSNRdB','Status','LeftSNRdB','RightSNRdB'});
    gains=cell2table(gain_rows,'VariableNames',{'Group','Comparison','TargetBER','ReferenceSNRdB','AxisSNRdB','MeasuredGainDb','BootstrapMeanDb', ...
        'BootstrapMedianDb','CI95LowerDb','CI95UpperDb','ResolvedResamples','Resamples','ResolutionFraction','Support','Magnitude'});
    writetable(summary,fullfile(output,'r22dcb_ber_summary.csv'));
    writetable(summary(:,{'Group','Method','SNRdB','Lower95','Upper95'}),fullfile(output,'r22dcb_ber_confidence_intervals.csv'));
    writetable(crossings,fullfile(output,'r22dcb_ber_snr_crossings.csv'));
    writetable(gains,fullfile(output,'r22dcb_ber_snr_gain_bootstrap.csv'));
    writetable(gains(startsWith(gains.Group,"VELOCITY_"),:),fullfile(output,'r22dcb_ber_velocity_summary.csv'));
    writetable(gains(gains.Group=="HIGH_MOBILITY",:),fullfile(output,'r22dcb_ber_high_mobility_summary.csv'));
    result=struct('Summary',summary,'Crossings',crossings,'Gains',gains);
end