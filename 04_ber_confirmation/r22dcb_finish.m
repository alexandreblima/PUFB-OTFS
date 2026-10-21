function result=r22dcb_finish(methods,p,output)
    if ischar(methods)&&strcmp(methods,'SELFTEST'),result=synthetic_check(p,output);return;end
    if ischar(methods)&&strcmp(methods,'A')
        loaded=load(fullfile(output,'r22dcb_frozen_methods.mat'),'methods');methods=loaded.methods;
        result=phase_a_statistics(methods,p,output);write_report(methods,p,output,result,[]);return
    end
    decision=jsondecode(fileread(fullfile(output,'r22dcb_phase_a_decision.json')));
    assert(strcmp(decision.R22D_CONFIRMATION,'CONFIRMED'));
    per_channel=r22dcb_math('read',fullfile(output,'r22dcb_ber_per_channel.csv'));
    gains=r22dcb_gains(methods,p,per_channel,output);
    full=table();
    for velocity=1:5
        loaded=load(fullfile(output,sprintf('full_control_v%d.mat',velocity)),'full_table');full=[full;loaded.full_table]; %#ok<AGROW>
    end
    shared=sortrows(full(string(full.Method)==methods(2).Name,:),{'VelocityIndex','ChannelID','SNRdB'});
    axis=sortrows(full(string(full.Method)==methods(3).Name,:),{'VelocityIndex','ChannelID','SNRdB'});
    assert(height(shared)==150&&height(axis)==150);
    difference=mean(reshape(axis.BER-shared.BER,3,50),1).';
    stats=r22dcb_math('bootstrap',difference,10,5,p.seeds.BER_BOOT,p.bootstrap_replicates);
    full.FullControlMeanDifference=repmat(stats.Mean,height(full),1);
    full.FullControlCI95Lower=repmat(stats.Lower,height(full),1);full.FullControlCI95Upper=repmat(stats.Upper,height(full),1);
    full.FullControlPass=repmat(stats.Lower<=0&&stats.Upper>=0,height(full),1);
    writetable(full,fullfile(output,'r22dcb_ber_full_receiver_control.csv'));
    gains.FullControlPass=stats.Lower<=0&&stats.Upper>=0;gains.FullControlStatistics=stats;
    result=write_report(methods,p,output,decision,gains);ber_figures(methods,p,gains,output);
end

function decision=phase_a_statistics(methods,p,output)
    mse=r22dcb_math('read',fullfile(output,'r22dcb_phase_a_exact_mse.csv'));
    full=r22dcb_math('read',fullfile(output,'r22dcb_phase_a_full_mse.csv'));
    delta=zeros(1000,8,4);full_values=delta;
    for method=1:4
        T=sortrows(mse(mse.Method==methods(method).Name,:),{'VelocityIndex','ChannelID','SNRdB'});
        assert(height(T)==8000&&all(T.RetainedEntries==p.K)&&isequal(T.VelocityIndex,repelem((1:5)',1600)) ...
            &&isequal(T.ChannelID,repmat(repelem((1:200)',8),5,1))&&isequal(T.SNRdB,repmat((3:3:24)',1000,1)));
        delta(:,:,method)=reshape(T.DeltaMSE,8,1000).';
        T=sortrows(full(full.Method==methods(method).Name,:),{'VelocityIndex','ChannelID','SNRdB'});
        full_values(:,:,method)=reshape(T.FullMSE,8,1000).';
    end
    invariance=max(abs(full_values-full_values(:,:,1)),[],'all');assert(invariance<=p.full_mse_tolerance);
    groups=struct('Name',{'AGGREGATE','LOW','MID','HIGH'},'Velocities',{1:5,1:5,1:5,1:5},'Points',{[],1:3,4:6,7:8});
    for velocity=1:5,groups(end+1)=struct('Name',sprintf('VELOCITY_%d',velocity),'Velocities',velocity,'Points',[]);end %#ok<AGROW>
    groups(end+1)=struct('Name','MACH5_21_24','Velocities',5,'Points',7:8);
    for point=1:8,groups(end+1)=struct('Name',sprintf('SNR_%d',3*point),'Velocities',1:5,'Points',point);end %#ok<AGROW>
    rows=cell(54,12);cursor=0;
    for reference=[2,1,4]
        for group=1:numel(groups)
            indices=reshape((1:1000)',200,5);indices=indices(:,groups(group).Velocities);indices=indices(:);
            left=r22dcb_math('channelmean',delta(:,:,reference),groups(group).Points);
            right=r22dcb_math('channelmean',delta(:,:,3),groups(group).Points);difference=right(indices)-left(indices);
            stats=r22dcb_math('bootstrap',difference,200,numel(groups(group).Velocities),p.seeds.MSE_BOOT,p.bootstrap_replicates);
            cursor=cursor+1;rows(cursor,:)={"AXIS_MINUS_"+methods(reference).Name,string(groups(group).Name),stats.Mean,stats.Lower,stats.Upper, ...
                mean(left(indices)),mean(right(indices)),(mean(left(indices))-mean(right(indices)))/mean(left(indices)), ...
                stats.Median,stats.SD,stats.Standardized,numel(indices)};
        end
    end
    paired=cell2table(rows,'VariableNames',{'Comparison','Group','MeanDifference','CI95Lower','CI95Upper','ReferenceMean','AxisMean', ...
        'RelativeReduction','PairedMedian','PairedSD','StandardizedEffect','PhysicalChannels'});
    writetable(paired,fullfile(output,'r22dcb_phase_a_bootstrap.csv'));
    writetable(paired(startsWith(paired.Group,"SNR_"),:),fullfile(output,'r22dcb_phase_a_snr_breakdown.csv'));
    writetable(paired(startsWith(paired.Group,"VELOCITY_"),:),fullfile(output,'r22dcb_phase_a_velocity_breakdown.csv'));
    writetable(paired(paired.Group=="MACH5_21_24",:),fullfile(output,'r22dcb_phase_a_stress_cell.csv'));
    primary=paired(paired.Comparison=="AXIS_MINUS_"+methods(2).Name&paired.Group=="AGGREGATE",:);
    confirmation='INCONCLUSIVE';if primary.CI95Upper<0,confirmation='CONFIRMED';elseif primary.MeanDifference>=0,confirmation='NOT_CONFIRMED';end
    decision=struct('R22D_CONFIRMATION',confirmation,'PHASE_A_CONFIRMATION',confirmation,'SharedDeltaMSE',primary.ReferenceMean, ...
        'AxisDeltaMSE',primary.AxisMean,'MeanDifference',primary.MeanDifference,'RelativeReduction',primary.RelativeReduction, ...
        'CI95Lower',primary.CI95Lower,'CI95Upper',primary.CI95Upper,'PairedMedian',primary.PairedMedian,'PairedSD',primary.PairedSD, ...
        'StandardizedPairedEffect',primary.StandardizedEffect,'FullReceiverInvarianceResidual',invariance,'PHASE_B_AUTHORIZED',false);
    r22a_io('json',fullfile(output,'r22dcb_phase_a_decision.json'),decision);
    fig=new_figure();tiledlayout(2,3,'TileSpacing','compact');labels={'D4','Shared','Independent'};
    for velocity=1:5
        nexttile;hold on;indices=(velocity-1)*200+(1:200);
        for method=1:3,semilogy(p.phase_a_snr_db,mean(delta(indices,:,method),1),'-o','DisplayName',labels{method});end
        grid on;set(gca,'YScale','log');title(p.velocity_labels{velocity});xlabel('SNR (dB)');ylabel('Exact excess MSE');
        if velocity==1,legend('Location','best');end
    end
    export_figure(fig,output,'fig_r22dcb_confirm_mse_vs_snr');
    selected=paired(paired.Comparison=="AXIS_MINUS_"+methods(2).Name&startsWith(paired.Group,"SNR_"),:);
    fig=new_figure();errorbar(p.phase_a_snr_db,selected.MeanDifference,selected.MeanDifference-selected.CI95Lower,selected.CI95Upper-selected.MeanDifference,'-o');
    hold on;yline(0);xlabel('SNR (dB)');ylabel('Independent - shared exact excess MSE, CI95');grid on;
    export_figure(fig,output,'fig_r22dcb_confirm_difference_vs_snr');
end

function flags=write_report(methods,p,output,decision,ber)
    authorized=~isempty(ber);gain1=NaN;gain_d4=NaN;gain4=NaN;support1='NOT_RESOLVED';support4='NOT_RESOLVED';high_support='NOT_RESOLVED';full_pass=false;
    interpretation='R22D_AXIS_SPECIFIC_GAIN_FAILED_STRICT_REPLICATION';
    if strcmp(decision.R22D_CONFIRMATION,'CONFIRMED'),interpretation='AXIS_SPECIFIC_MSE_GAIN_CONFIRMED_WITHOUT_RESOLVED_BER_SNR_GAIN';end
    if authorized
        primary=gain_row(ber.Gains,'AGGREGATE',methods(2).Name,1e-3);fixed=gain_row(ber.Gains,'AGGREGATE',methods(1).Name,1e-3);
        secondary=gain_row(ber.Gains,'AGGREGATE',methods(2).Name,1e-4);high=gain_row(ber.Gains,'HIGH_MOBILITY',methods(2).Name,1e-3);
        low=gain_row(ber.Gains,'LOW_MOBILITY',methods(2).Name,1e-3);
        gain1=primary.MeasuredGainDb;gain_d4=fixed.MeasuredGainDb;gain4=secondary.MeasuredGainDb;
        support1=char(primary.Support);support4=char(secondary.Support);high_support=char(high.Support);full_pass=ber.FullControlPass;
        if strcmp(support1,'SUPPORTED')
            interpretation='AXIS_SPECIFIC_RECEIVER_AWARE_PUFB_CONFIRMED_WITH_OPERATIONAL_BER_SNR_GAIN';
            if strcmp(high_support,'SUPPORTED')&&low.Support~="SUPPORTED"&&isfinite(low.MeasuredGainDb)&&high.MeasuredGainDb>low.MeasuredGainDb
                interpretation='AXIS_SPECIFIC_OPERATIONAL_GAIN_IS_HIGH_MOBILITY_DOMINANT';
            end
        end
    end
    flags=struct('ALPHA_NU',methods(3).AlphaNu,'ALPHA_TAU',methods(3).AlphaTau,'ALPHA_SHARED',methods(2).AlphaNu, ...
        'PRIMARY_RETAINED_ENTRIES',p.K,'NEW_TRANSFORM_TRAINED',false,'SUPPORT_REBUILT',false,'FILTER_LENGTH_ESCALATED',false, ...
        'TEST_USED_FOR_SELECTION',false,'PHASE_A_CHANNELS',1000,'PHASE_A_CONFIRMATION',decision.R22D_CONFIRMATION, ...
        'PHASE_B_AUTHORIZED',authorized,'BER_CHANNELS',500*authorized,'BER_SNR_GRID',0:24,'BER_TARGET_PRIMARY',1e-3, ...
        'BER_TARGET_SECONDARY',1e-4,'BER_SNR_GAIN_1E3',support1,'BER_SNR_GAIN_1E4',support4, ...
        'HIGH_MOBILITY_BER_GAIN',high_support,'DELTA_SNR_AXIS_VS_SHARED_1E3_DB',gain1,'DELTA_SNR_AXIS_VS_D4_1E3_DB',gain_d4, ...
        'DELTA_SNR_AXIS_VS_SHARED_1E4_DB',gain4,'FULL_RECEIVER_CONTROL_PASS',full_pass, ...
        'MANUSCRIPT_MODIFIED',false,'PRIMARY_INTERPRETATION',interpretation);
    r22a_io('json',fullfile(output,'r22dcb_decision_flags.json'),flags);
    paired=r22dcb_math('read',fullfile(output,'r22dcb_phase_a_bootstrap.csv'));
    lines=["# R2.2D-CONFIRM + BER";"## Phase A"; ...
        sprintf('1. alpha_nu: %.15g.',methods(3).AlphaNu);sprintf('2. alpha_tau: %.15g.',methods(3).AlphaTau); ...
        sprintf('3. alpha_shared: %.15g.',methods(2).AlphaNu);"4. K: 33178.";"5. Confirmation channels: 1000 (200 per velocity)."; ...
        sprintf('6. Shared DeltaMSE: %.12g.',decision.SharedDeltaMSE);sprintf('7. Axis DeltaMSE: %.12g.',decision.AxisDeltaMSE); ...
        sprintf('8. Axis minus shared: %.12g.',decision.MeanDifference);sprintf('9. Relative reduction: %.9g%%.',100*decision.RelativeReduction); ...
        sprintf('10. CI95: [%.12g,%.12g].',decision.CI95Lower,decision.CI95Upper); ...
        sprintf('11. Standardized paired effect %.12g; median %.12g; SD %.12g.',decision.StandardizedPairedEffect,decision.PairedMedian,decision.PairedSD)];
    labels={'LOW','MID','HIGH','VELOCITY_1','VELOCITY_2','VELOCITY_3','VELOCITY_4','VELOCITY_5','MACH5_21_24'};
    for index=1:numel(labels)
        selected=paired(paired.Comparison=="AXIS_MINUS_"+methods(2).Name&paired.Group==string(labels{index}),:);
        lines(end+1)=sprintf('%d. %s difference %.12g CI95 [%.12g,%.12g].',11+index,labels{index},selected.MeanDifference,selected.CI95Lower,selected.CI95Upper); %#ok<AGROW>
    end
    lines=[lines;sprintf('21. Full receiver invariance residual: %.12g.',decision.FullReceiverInvarianceResidual); ...
        "22. PHASE_A_CONFIRMATION="+string(decision.R22D_CONFIRMATION)];
    if authorized
        counts=r22dcb_math('read',fullfile(output,'r22dcb_ber_counts.csv'));
        unique_conditions=counts(string(counts.Method)==methods(1).Name,:);total_bits=sum(unique_conditions.Bits);
        lines=[lines;"## Phase B";"23. BER channels: 500 (100 per velocity), independent of Phase A."; ...
            sprintf('24. Total paired transmitted bits: %.0f per method; %.0f counted across four methods.',total_bits,4*total_bits)];
        for method=1:3
            T=ber.Summary(ber.Summary.Group=="AGGREGATE"&ber.Summary.Method==methods(method).Name,:);
            lines(end+1)=sprintf('%d. %s aggregate curve:25 measured points,0:24dB; min measured BER %.12g; no extrapolation.',24+method,methods(method).Name,min(T.BER)); %#ok<AGROW>
        end
        low_target=gain_row(ber.Gains,'AGGREGATE',methods(2).Name,1e-2);
        lines=[lines;sprintf('28. Shared SNR at1e-2: %.12g dB.',low_target.ReferenceSNRdB); ...
            sprintf('29. Axis SNR at1e-2: %.12g dB.',low_target.AxisSNRdB);sprintf('30. DeltaSNR at1e-2: %.12g dB (%s).',low_target.MeasuredGainDb,low_target.Support); ...
            sprintf('31. Shared SNR at1e-3: %.12g dB.',primary.ReferenceSNRdB);sprintf('32. Axis SNR at1e-3: %.12g dB.',primary.AxisSNRdB); ...
            sprintf('33. DeltaSNR axis-vs-shared at1e-3: %.12g dB (%s).',primary.MeasuredGainDb,primary.Support); ...
            sprintf('34. CI95 DeltaSNR1e-3 [%.12g,%.12g]; resolved bootstrap %d/%d.',primary.CI95LowerDb,primary.CI95UpperDb,primary.ResolvedResamples,primary.Resamples); ...
            sprintf('35. DeltaSNR axis-vs-D4 at1e-3: %.12g dB (%s).',fixed.MeasuredGainDb,fixed.Support); ...
            sprintf('36. Shared SNR at1e-4: %.12g dB.',secondary.ReferenceSNRdB);sprintf('37. Axis SNR at1e-4: %.12g dB.',secondary.AxisSNRdB); ...
            sprintf('38. DeltaSNR1e-4: %.12g dB (%s).',secondary.MeasuredGainDb,secondary.Support)];
        for velocity=1:5
            selected=gain_row(ber.Gains,sprintf('VELOCITY_%d',velocity),methods(2).Name,1e-3);
            lines(end+1)=sprintf('%d. %s DeltaSNR1e-3 %.12g dB CI95 [%.12g,%.12g] (%s).',38+velocity,p.velocity_labels{velocity}, ...
                selected.MeasuredGainDb,selected.CI95LowerDb,selected.CI95UpperDb,selected.Support); %#ok<AGROW>
        end
        lines=[lines;sprintf('44. High mobility DeltaSNR1e-3 %.12g dB CI95 [%.12g,%.12g].',high.MeasuredGainDb,high.CI95LowerDb,high.CI95UpperDb); ...
            "45. BER_SNR_GAIN_1E3="+string(support1);"46. BER_SNR_GAIN_1E4="+string(support4);"47. HIGH_MOBILITY_BER_GAIN="+string(high_support); ...
            "48. FULL_RECEIVER_CONTROL_PASS="+string(full_pass);"49. PRIMARY_INTERPRETATION="+string(interpretation); ...
            "50. NEW_TRANSFORM_TRAINED=FALSE";"51. SUPPORT_REBUILT=FALSE";"52. FILTER_LENGTH_ESCALATED=FALSE"; ...
            "53. TEST_USED_FOR_SELECTION=FALSE";"54. MANUSCRIPT_MODIFIED=FALSE"; ...
            sprintf('Full-channel BER diagnostic difference %.12g CI95 [%.12g,%.12g].',ber.FullControlStatistics.Mean,ber.FullControlStatistics.Lower,ber.FullControlStatistics.Upper)];
    else
        lines(end+1)="Phase B not executed. Authorization requires persisted CONFIRMED Phase-A decision.";
    end
    lines=[lines;"";"## Statistical And Procedural Notes"; ...
        "All transforms and exact R2.2D TRAIN masks loaded, never rebuilt. Pre-generation receipt freezes BOTH phases, sources, seeds, stopping, targets, bootstrap and report code. Phase-B authorization is persisted before its independent channel bank."; ...
        "Phase A equally weights LOW/MID/HIGH and velocities. Phase-B curves equally weight per-velocity count ratios. Physical-channel bootstrap preserves paired methods and all SNRs/frames within a channel;10000 resamples."; ...
        "Crossings use only two adjacent positive measured BER points with at least300 aggregate errors each. No polynomial/spline/extrapolation. Unresolved targets remain NaN in CSV/report (JSON null), never assigned zero gain. At least95% resamples must resolve; reported percentiles are conditional on resolved resamples with counts disclosed."; ...
        "Measured aggregate Wilson intervals are averages of Bonferroni simultaneous velocity intervals; method comparisons use paired channel bootstrap, not binomial independence of fading frames."; ...
        "All four methods share rounds and data/noise; stop uses all three primary methods>=300errors OR hard cap. Batch rounding may overshoot300errors but never the5e7bits cap. Secondary method never changes stopping."; ...
        "Full trace-MSE unitary invariance does not imply full-channel QPSK BER invariance: different unitary precoders change finite-constellation symbol geometry. Full BER control is empirical; no forced equality and no selection based on its result."; ...
        "MSE confirmation does not by itself establish operational BER gain; no unresolved target is extrapolated. No alpha/support/K/receiver/channel/grid changes or manuscript edits."];
    r22a_io('text',fullfile(output,'r22dcb_final_report.md'),join(lines,newline));
    r22a_io('text',fullfile(output,'r22dcb_executive_summary.md'),join(lines,newline));
end

function row=gain_row(T,group,reference,target)
    row=T(T.Group==string(group)&T.Comparison=="AXIS_VS_"+string(reference)&T.TargetBER==target,:);assert(height(row)==1);
end

function ber_figures(methods,~,ber,output)
    plot_curve(methods,ber.Summary,'AGGREGATE',output,'fig_r22dcb_ber_aggregate');
    plot_curve(methods,ber.Summary,'HIGH_MOBILITY',output,'fig_r22dcb_ber_high_mobility');
    plot_curve(methods,ber.Summary,'VELOCITY_5',output,'fig_r22dcb_ber_mach5');
    T=ber.Gains(ber.Gains.Comparison=="AXIS_VS_"+methods(2).Name&ber.Gains.TargetBER==1e-3&startsWith(ber.Gains.Group,"VELOCITY_"),:);
    fig=new_figure();resolved=isfinite(T.MeasuredGainDb)&T.Support~="NOT_RESOLVED";
    if any(resolved)
        positions=find(resolved);errorbar(positions,T.MeasuredGainDb(resolved),T.MeasuredGainDb(resolved)-T.CI95LowerDb(resolved), ...
            T.CI95UpperDb(resolved)-T.MeasuredGainDb(resolved),'o');hold on;yline(0);
    end
    if any(~resolved),text(find(~resolved),zeros(nnz(~resolved),1),repmat({'NOT RESOLVED'},nnz(~resolved),1),'Rotation',45);end
    xlim([.5,5.5]);xticks(1:5);xticklabels({'50 km/h','100 km/h','500 km/h','Mach 1','Mach 5'});ylabel('Independent vs. shared SNR gain at BER1e-3 (dB)');grid on;
    export_figure(fig,output,'fig_r22dcb_snr_gain_1e3_by_velocity');
    fig=new_figure();tiledlayout(1,3,'TileSpacing','compact');targets=[1e-2,1e-3,1e-4];
    for target=targets
        nexttile;T=ber.Gains(ber.Gains.Group=="AGGREGATE"&ber.Gains.TargetBER==target& ...
            ismember(ber.Gains.Comparison,["AXIS_VS_"+methods(2).Name,"AXIS_VS_"+methods(1).Name]),:);
        resolved=isfinite(T.MeasuredGainDb)&T.Support~="NOT_RESOLVED";
        if any(resolved)
            positions=find(resolved);errorbar(positions,T.MeasuredGainDb(resolved),T.MeasuredGainDb(resolved)-T.CI95LowerDb(resolved), ...
                T.CI95UpperDb(resolved)-T.MeasuredGainDb(resolved),'o');hold on;yline(0);
        end
        if any(~resolved),text(find(~resolved),zeros(nnz(~resolved),1),repmat({'NOT RESOLVED'},nnz(~resolved),1),'Rotation',45);end
        xlim([.5,2.5]);xticks(1:2);xticklabels({'vs shared','vs D4'});ylabel('SNR gain (dB)');title(sprintf('BER %.0e',target));grid on;
    end
    export_figure(fig,output,'fig_r22dcb_snr_gain_summary');
end

function plot_curve(methods,summary,group,output,name)
    fig=new_figure();hold on;labels={'D4','Shared','Independent','Compressibility-trained shared'};
    for method=1:4
        T=summary(summary.Group==string(group)&summary.Method==methods(method).Name,:);positive=T.BER>0;
        plotted=T.BER;plotted(~positive)=NaN;
        errorbar(T.SNRdB,plotted,max(0,T.BER-T.Lower95),max(0,T.Upper95-T.BER),'-o','DisplayName',labels{method});
        if any(~positive),scatter(T.SNRdB(~positive),T.Upper95(~positive),30,'v','DisplayName',[labels{method},' zero-error upper bound']);end
    end
    set(gca,'YScale','log');grid on;xlabel('SNR (dB)');ylabel('Measured BER, confidence interval');legend('Location','best');title(strrep(group,'_',' '));
    export_figure(fig,output,name);
end

function fig=new_figure()
    fig=figure('Visible','off','Color','w','Position',[100,100,1100,650]);
end

function export_figure(fig,output,name)
    exportgraphics(fig,fullfile(output,[name,'.pdf']),'ContentType','vector');
    exportgraphics(fig,fullfile(output,[name,'.png']),'Resolution',180);close(fig);
end

function result=synthetic_check(p,output)
    smoke=fullfile(output,'synthetic_preflight');if ~isfolder(smoke),mkdir(smoke);end
    loaded=load(fullfile(output,'r22dcb_frozen_methods.mat'),'methods');methods=loaded.methods;
    save(fullfile(smoke,'r22dcb_frozen_methods.mat'),'methods');
    rows=cell(32000,8);full_rows=cell(32000,5);cursor=0;
    for velocity=1:5
        for channel=1:200
            for method=1:4
                for point=1:8
                    cursor=cursor+1;values=[.05,.03,.02,.035];delta=values(method)+channel*1e-7;
                    rows(cursor,:)={methods(method).Name,velocity,channel,point*3,p.K,.1+delta,.1,delta};
                    full_rows(cursor,:)={methods(method).Name,velocity,channel,point*3,.1};
                end
            end
        end
    end
    mse=cell2table(rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','RetainedEntries','MSESparse','MSEFull','DeltaMSE'});
    full=cell2table(full_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','FullMSE'});
    writetable(mse,fullfile(smoke,'r22dcb_phase_a_exact_mse.csv'));writetable(full,fullfile(smoke,'r22dcb_phase_a_full_mse.csv'));
    decision=r22dcb_finish('A',p,smoke);assert(strcmp(decision.R22D_CONFIRMATION,'CONFIRMED'));
    rows=cell(50000,9);count_rows=cell(500,13);cursor=0;count_cursor=0;shifts=[-1,0,2,-.2];
    for velocity=1:5
        for snr=0:24
            for method=1:4
                total=0;
                for channel=1:100
                    rate=.1*10^(-(snr+shifts(method))/6)*(1+.001*channel);errors=round(1000000*rate);total=total+errors;cursor=cursor+1;
                    rows(cursor,:)={methods(method).Name,velocity,channel,snr,errors,1000000,1736,p.K,p.seeds.B+1000*velocity+channel};
                end
                interval=r22dcb_math('wilson',total,100000000);count_cursor=count_cursor+1;
                count_rows(count_cursor,:)={methods(method).Name,velocity,snr,total,100000000,1736,100,false,true,total/100000000,interval(1),interval(2),p.K};
            end
        end
        control_rows=cell(60,8);control_cursor=0;
        for snr=[6,15,24]
            for channel=1:10
                for method=[2,3]
                    control_cursor=control_cursor+1;control_rows(control_cursor,:)={methods(method).Name,velocity,channel,snr,100,100000,1736,.001};
                end
            end
        end
        full_table=cell2table(control_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','Errors','Bits','Frames','BER'});
        save(fullfile(smoke,sprintf('full_control_v%d.mat',velocity)),'full_table');
    end
    channels=cell2table(rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','Errors','Bits','Frames','RetainedEntries','Seed'});
    counts=cell2table(count_rows,'VariableNames',{'Method','VelocityIndex','SNRdB','Errors','Bits','FramesPerChannel','ChannelsUsed','ZeroError','CapReached','BER','WilsonLower95','WilsonUpper95','RetainedEntries'});
    writetable(channels,fullfile(smoke,'r22dcb_ber_per_channel.csv'));writetable(counts,fullfile(smoke,'r22dcb_ber_counts.csv'));
    flags=r22dcb_finish(methods,p,smoke);assert(strcmp(flags.BER_SNR_GAIN_1E3,'SUPPORTED')&&flags.FULL_RECEIVER_CONTROL_PASS ...
        &&abs(flags.DELTA_SNR_AXIS_VS_SHARED_1E3_DB-2)<.01);
    result=struct('PASS',true,'PHYSICAL_CHANNELS_GENERATED',false,'PHASE_A_SYNTHETIC_CONFIRMED',true,'KNOWN_2DB_GAIN_REPRODUCED',true);
    r22a_io('json',fullfile(output,'r22dcb_preflight_report.json'),result);
end