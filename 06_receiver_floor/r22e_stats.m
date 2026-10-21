function result=r22e_stats(action,p,output)
    if strcmp(action,'MSE'),result=mse_stats(p,output);else,result=ber_stats(p,output);end
end

function result=mse_stats(p,output)
    mse=r22e_math('read',fullfile(output,'r22e_test_exact_mse.csv'));diagnostic=r22e_math('read',fullfile(output,'r22e_test_mismatch_decomposition.csv'));
    names=unique(mse.Variant,'stable');matrices=struct();
    for index=1:numel(names)
        T=sortrows(mse(mse.Variant==names(index),:),{'VelocityIndex','ChannelID','SNRdB'});
        assert(height(T)==17000&&isequal(T.SNRdB,repmat((20:36)',1000,1))&&all(T.RetainedEntries==33178));
        matrices.(char(names(index)))=reshape(T.MSE,17,1000).';
    end
    pairs={'AXIS_REG','AXIS_ETA0';'AXIS_REG','SHARED_REG';'AXIS_COMMON','SHARED_COMMON'};
    regions={20:36,20:24,25:29,30:36};labels={'AGGREGATE','LOW','MID','HIGH'};rows=cell(12,9);cursor=0;
    for pair=1:3
        for region=1:4
            selected=ismember(p.development_snr,regions{region});left=mean(matrices.(pairs{pair,2})(:,selected),2);right=mean(matrices.(pairs{pair,1})(:,selected),2);
            stats=r22e_math('bootstrap',right-left,200,p.seeds.BOOTSTRAP);cursor=cursor+1;
            rows(cursor,:)={string(pairs{pair,1})+"_MINUS_"+string(pairs{pair,2}),string(labels{region}),mean(left),mean(right),stats.Mean,stats.Lower,stats.Upper,stats.Median,stats.Standardized};
        end
    end
    comparisons=cell2table(rows,'VariableNames',{'Comparison','Region','ReferenceMean','CandidateMean','MeanDifference','CI95Lower','CI95Upper','PairedMedian','StandardizedEffect'});
    writetable(comparisons,fullfile(output,'r22e_test_mse_bootstrap.csv'));
    metrics={'MSE','T_REG','T_OMIT','T_CROSS','T_NOISE','W_F','WSE_F','CondA','MismatchPenalty'};rows=cell(9,6);
    for index=1:9
        field=metrics{index};base=sortrows(diagnostic(diagnostic.Variant=="AXIS_ETA0",:),{'VelocityIndex','ChannelID','SNRdB'});
        reg=sortrows(diagnostic(diagnostic.Variant=="AXIS_REG",:),{'VelocityIndex','ChannelID','SNRdB'});
        left=reshape(base.(field),17,1000).';right=reshape(reg.(field),17,1000).';selected=p.development_snr>=30;
        left=mean(left(:,selected),2);right=mean(right(:,selected),2);stats=r22e_math('bootstrap',right-left,200,p.seeds.BOOTSTRAP);
        rows(index,:)={string(field),mean(left),mean(right),stats.Mean,stats.Lower,stats.Upper};
    end
    mechanism=cell2table(rows,'VariableNames',{'Metric','ReferenceMean','CandidateMean','MeanDifference','CI95Lower','CI95Upper'});
    writetable(mechanism,fullfile(output,'r22e_mechanism_bootstrap.csv'));
    result=struct('Comparisons',comparisons,'Mechanism',mechanism,'Matrices',matrices);
end

function result=ber_stats(p,output)
    channels=r22e_math('read',fullfile(output,'r22e_ber_per_channel.csv'));names=["AXIS_REG","SHARED_REG","D4_REG","R22A_REG","AXIS_ETA0"];
    errors=zeros(200,37,5,5);bits=errors;
    for method=1:5
        for velocity=1:5
            T=sortrows(channels(channels.Variant==names(method)&channels.VelocityIndex==velocity,:),{'ChannelID','SNRdB'});
            assert(height(T)==7400&&isequal(T.ChannelID,repelem((1:200)',37))&&isequal(T.SNRdB,repmat((0:36)',200,1)));
            errors(:,:,method,velocity)=reshape(T.Errors,37,200).';bits(:,:,method,velocity)=reshape(T.Bits,37,200).';
        end
    end
    groups=struct('Name',{'AGGREGATE','HIGH_MOBILITY'},'Velocities',{1:5,[4,5]});
    for velocity=1:5,groups(end+1)=struct('Name',sprintf('VELOCITY_%d',velocity),'Velocities',velocity);end %#ok<AGROW>
    summary_rows=cell(7*5*37,9);cross_rows=cell(7*5*3,7);gain_rows=cell(7*3*3,14);shape_rows=cell(7*5,8);high_rows=cell(7,4);
    summary_cursor=0;cross_cursor=0;gain_cursor=0;shape_cursor=0;
    for group=1:7
        velocities=groups(group).Velocities;rates=zeros(5,37);lower=rates;upper=rates;totals=rates;total_bits=rates;
        boot=zeros(10000,37,5);stream=RandStream('mt19937ar','Seed',p.seeds.BOOTSTRAP);
        for velocity=velocities
            draws=randi(stream,200,[200,10000]);weights=accumarray([draws(:),repelem((1:10000)',200)],1,[200,10000]);
            for method=1:5
                current_errors=errors(:,:,method,velocity);current_bits=bits(:,:,method,velocity);
                totals(method,:)=totals(method,:)+sum(current_errors,1);total_bits(method,:)=total_bits(method,:)+sum(current_bits,1);
                rates(method,:)=rates(method,:)+sum(current_errors,1)./sum(current_bits,1)/numel(velocities);
                boot(:,:,method)=boot(:,:,method)+(weights'*current_errors)./(weights'*current_bits)/numel(velocities);
                for point=1:37
                    interval=r22dcb_math('wilson',sum(current_errors(:,point)),sum(current_bits(:,point)),.05/numel(velocities));
                    lower(method,point)=lower(method,point)+interval(1)/numel(velocities);upper(method,point)=upper(method,point)+interval(2)/numel(velocities);
                end
            end
        end
        eligible=totals>=500;crossings=nan(5,3);sampled_crossings=nan(10000,5,3);
        for method=1:5
            for point=1:37
                summary_cursor=summary_cursor+1;summary_rows(summary_cursor,:)={string(groups(group).Name),names(method),point-1,rates(method,point), ...
                    lower(method,point),upper(method,point),totals(method,point),totals(method,point)/total_bits(method,point),eligible(method,point)};
            end
            for target_index=1:3
                cross=r22dcb_math('crossing',0:36,rates(method,:),p.ber_targets(target_index),eligible(method,:));
                crossings(method,target_index)=cross.SNR;cross_cursor=cross_cursor+1;
                cross_rows(cross_cursor,:)={string(groups(group).Name),names(method),p.ber_targets(target_index),cross.SNR,string(cross.Status),cross.LeftIndex-1,cross.RightIndex-1};
                if isfinite(cross.SNR)
                    for replicate=1:10000
                        sampled=r22dcb_math('crossing',0:36,boot(replicate,:,method),p.ber_targets(target_index),eligible(method,:));
                        sampled_crossings(replicate,method,target_index)=sampled.SNR;
                    end
                end
            end
            selected=(0:36)>=30&rates(method,:)>0;slope=NaN;shape='UNRESOLVED';
            if nnz(selected)>=3
                snrs=0:36;fitted=polyfit(snrs(selected),log10(rates(method,selected)),1);slope=fitted(1);
                if slope<=-.03,shape='CLEAR_DESCENT';elseif slope>=.03,shape='UPTURN';else,shape='APPROX_FLAT';end
            end
            [minimum,min_index]=min(rates(method,:));shape_cursor=shape_cursor+1;
            shape_rows(shape_cursor,:)={string(groups(group).Name),names(method),minimum,min_index-1,rates(method,31),rates(method,34),rates(method,37),slope};
            shape_classes(shape_cursor,1)=string(shape); %#ok<AGROW>
        end
        for reference=[2,5,3]
            for target_index=1:3
                difference=sampled_crossings(:,reference,target_index)-sampled_crossings(:,1,target_index);valid=isfinite(difference);
                measured=crossings(reference,target_index)-crossings(1,target_index);lo=NaN;hi=NaN;average=NaN;middle=NaN;support='NOT_RESOLVED';
                if isfinite(measured)&&mean(valid)>=.95
                    samples=sort(difference(valid));lo=samples(ceil(.025*numel(samples)));hi=samples(ceil(.975*numel(samples)));
                    average=mean(samples);middle=median(samples);support='NOT_SUPPORTED';if lo>0,support='SUPPORTED';end
                end
                gain_cursor=gain_cursor+1;gain_rows(gain_cursor,:)={string(groups(group).Name),"AXIS_REG_VS_"+names(reference),p.ber_targets(target_index), ...
                    crossings(reference,target_index),crossings(1,target_index),measured,average,middle,lo,hi,nnz(valid),10000,mean(valid),string(support)};
            end
        end
        difference=mean(boot(:,31:37,1)-boot(:,31:37,5),2);samples=sort(difference);
        high_rows(group,:)={string(groups(group).Name),mean(rates(1,31:37)-rates(5,31:37)),samples(250),samples(9750)};
        fprintf('R22E BER paired crossings group%s done\n',groups(group).Name);
    end
    summary=cell2table(summary_rows,'VariableNames',{'Group','Variant','SNRdB','BER','Lower95','Upper95','Errors','RawPooledBER','CrossingEligible'});
    crossings=cell2table(cross_rows,'VariableNames',{'Group','Variant','TargetBER','CrossingSNRdB','Status','LeftSNRdB','RightSNRdB'});
    gains=cell2table(gain_rows,'VariableNames',{'Group','Comparison','TargetBER','ReferenceSNRdB','AxisSNRdB','MeasuredGainDb','BootstrapMeanDb', ...
        'BootstrapMedianDb','CI95LowerDb','CI95UpperDb','ResolvedResamples','Resamples','ResolutionFraction','Support'});
    shapes=cell2table(shape_rows,'VariableNames',{'Group','Variant','MinimumBER','MinimumBERSNR','BER30','BER33','BER36','HighSNRLogSlope'});shapes.Shape=shape_classes;
    high=cell2table(high_rows,'VariableNames',{'Group','AxisREGMinusETA0HighBER','CI95Lower','CI95Upper'});
    writetable(summary,fullfile(output,'r22e_ber_summary.csv'));writetable(summary(:,{'Group','Variant','SNRdB','Lower95','Upper95'}),fullfile(output,'r22e_ber_confidence_intervals.csv'));
    writetable(crossings,fullfile(output,'r22e_ber_crossings.csv'));writetable(gains,fullfile(output,'r22e_snr_gain_bootstrap.csv'));
    writetable(shapes,fullfile(output,'r22e_highsnr_shape.csv'));writetable(gains(startsWith(gains.Group,"VELOCITY_"),:),fullfile(output,'r22e_velocity_summary.csv'));
    writetable(high,fullfile(output,'r22e_highsnr_ber_bootstrap.csv'));
    result=struct('Summary',summary,'Crossings',crossings,'Gains',gains,'Shapes',shapes,'HighBER',high);
end