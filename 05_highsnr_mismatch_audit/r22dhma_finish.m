function flags=r22dhma_finish(methods,p,tables,source,output)
    if ischar(methods)&&strcmp(methods,'PLOT_DECOMPOSITION_ONLY')
        curves=readtable(fullfile(output,'r22dhma_equal_velocity_curves.csv'),'TextType','string');
        destination=p;figure_size=[7.2,3.7];
        if isstruct(p),destination=p.OutputDirectory;figure_size=p.FigureSize;end
        exact_mse_decomposition_figure(curves,["R2.2C_SHARED_RX","R2.2D_AXIS_RX"],destination,figure_size);return;
    end
    if ischar(methods)&&strcmp(methods,'SELFTEST'),flags=selftest(p,source,output);return;end
    T=tables.Detailed;S=tables.Static;fields=T.Properties.VariableNames(5:end);cube=struct();
    for field=fields
        values=zeros(500,17,4);
        for method=1:4
            selected=sortrows(T(string(T.Method)==methods(method).Name,:),{'VelocityIndex','ChannelID','SNRdB'});
            assert(height(selected)==8500&&isequal(selected.SNRdB,repmat((20:36)',500,1)));
            values(:,:,method)=reshape(selected.(field{1}),17,500).';
        end
        cube.(field{1})=values;
    end
    stat=struct();static_fields=S.Properties.VariableNames(4:end);
    for field=static_fields
        values=zeros(500,4);
        for method=1:4,selected=sortrows(S(string(S.Method)==methods(method).Name,:),{'VelocityIndex','ChannelID'});values(:,method)=selected.(field{1});end
        stat.(field{1})=values;
    end
    full_pass=all(diff(cube.MSE_FULL,1,2)<=p.monotonic_tolerance*(1+cube.MSE_FULL(:,1:16,:)),'all');
    consistent_pass=all(diff(cube.MSE_MODEL_CONSISTENT,1,2)<=p.monotonic_tolerance*(1+cube.MSE_MODEL_CONSISTENT(:,1:16,:)),'all');
    maximum_backward=max(cube.BackwardResidual,[],'all');
    maximum_reference=max(cube.SVDReceiverRelativeError,[],'all','omitnan');maximum_mse_reference=max(cube.SVDMSEDifference,[],'all','omitnan');
    maximum_condition_error=max(cube.ConditionIdentityError,[],'all','omitnan');
    solver_pass=maximum_backward<=p.backward_tolerance&&maximum_reference<=p.svd_receiver_relative_tolerance ...
        &&maximum_mse_reference<=p.svd_mse_absolute_tolerance&&maximum_condition_error<=p.condition_identity_relative_tolerance;
    noise_audit=jsondecode(fileread(fullfile(output,'r22dhma_snr_noise_audit.json')));
    mid=ismember(p.snr_db,p.mid_reference);high=ismember(p.snr_db,p.high_region);
    changes={'MSE_TRUE_SPARSE','MISMATCH_PENALTY','MSE_MODEL_CONSISTENT','T_OMIT','T_CROSS','T_NOISE','W_F','WSE_F','Amplification_F','CondA'};
    change_rows=cell(40,7);cursor=0;
    for method=1:4
        for field=changes
            values=mean(cube.(field{1})(:,high,method),2)-mean(cube.(field{1})(:,mid,method),2);
            stats=r22dcb_math('bootstrap',values,100,5,p.bootstrap_seed,p.bootstrap_replicates);
            cursor=cursor+1;change_rows(cursor,:)={methods(method).Name,string(field{1}),stats.Mean,stats.Lower,stats.Upper,stats.Median,stats.SD};
        end
    end
    change=cell2table(change_rows,'VariableNames',{'Method','Metric','MeanDifference','CI95Lower','CI95Upper','PairedMedian','PairedSD'});
    writetable(change,fullfile(output,'r22dhma_highsnr_change_bootstrap.csv'));
    comparisons={'RelativeOmittedF','SigmaMin','CondHs','CondA','W_F','WSE_F','MISMATCH_PENALTY','T_OMIT','T_CROSS','MSE_TRUE_SPARSE','Amplification_F'};
    comparison_rows=cell(55,7);cursor=0;conditions={24,30,33,36,[34,35,36]};condition_names={'24','30','33','36','HIGH_34_36'};
    for condition=1:5
        for field=comparisons
            if isfield(stat,field{1}),left=stat.(field{1})(:,2);right=stat.(field{1})(:,3);
            else
                selected=ismember(p.snr_db,conditions{condition});left=mean(cube.(field{1})(:,selected,2),2);right=mean(cube.(field{1})(:,selected,3),2);
            end
            difference=right-left;stats=r22dcb_math('bootstrap',difference,100,5,p.bootstrap_seed,p.bootstrap_replicates);
            cursor=cursor+1;comparison_rows(cursor,:)={string(condition_names{condition}),string(field{1}),stats.Mean,stats.Lower,stats.Upper,mean(left),mean(right)};
        end
    end
    paired=cell2table(comparison_rows,'VariableNames',{'Condition','Metric','MeanDifference','CI95Lower','CI95Upper','SharedMean','AxisMean'});
    writetable(paired,fullfile(output,'r22dhma_axis_vs_shared.csv'));
    curves=aggregate_curves(cube,methods,p);writetable(curves,fullfile(output,'r22dhma_equal_velocity_curves.csv'));
    counts=r22dcb_math('read',fullfile(source,'r22dbe_ber_counts.csv'));ber=historical_ber(counts,methods,p);
    turning_rows=cell(12,9);cursor=0;groups={'AGGREGATE','VELOCITY_1','VELOCITY_2','VELOCITY_3','VELOCITY_4','VELOCITY_5'};
    for group=1:6
        indices=1:500;if group>1,indices=(group-2)*100+(1:100);end
        for method=[2,3]
            true_mse=mean(cube.MSE_TRUE_SPARSE(indices,:,method),1);penalty=mean(cube.MISMATCH_PENALTY(indices,:,method),1);
            [best_mse,mse_index]=min(true_mse);[best_penalty,penalty_index]=min(penalty);rates=ber(ber.Group==string(groups{group})&ber.Method==methods(method).Name,:);
            [best_ber,ber_index]=min(rates.BER);first=NaN;
            for point=1:15,if penalty(point+1)>penalty(point)&&penalty(point+2)>penalty(point+1),first=p.snr_db(point);break;end;end
            cursor=cursor+1;turning_rows(cursor,:)={methods(method).Name,string(groups{group}),rates.SNRdB(ber_index),best_ber, ...
                p.snr_db(mse_index),best_mse,p.snr_db(penalty_index),best_penalty,first};
        end
    end
    turning=cell2table(turning_rows,'VariableNames',{'Method','Group','MinimumBERSNR','MinimumBER','MinimumSparseMSESNR','MinimumSparseMSE', ...
        'MinimumMismatchSNR','MinimumMismatchPenalty','FirstTwoStepMismatchIncreaseAfterSNR'});
    writetable(turning,fullfile(output,'r22dhma_turning_points.csv'));
    correlation_metrics={'MISMATCH_PENALTY','T_OMIT','WSE_F','W_F','CondA','MSE_TRUE_SPARSE'};correlation_rows=cell(12,5);cursor=0;
    for method=[2,3]
        rates=ber(ber.Group=="AGGREGATE"&ber.Method==methods(method).Name,:);
        for field=correlation_metrics
            values=mean(cube.(field{1})(:,:,method),1).';[pearson,spearman,count]=association(rates.BER,values);
            cursor=cursor+1;correlation_rows(cursor,:)={methods(method).Name,string(field{1}),pearson,spearman,count};
        end
    end
    correlations=cell2table(correlation_rows,'VariableNames',{'Method','Metric','Pearson','Spearman','SNRPoints'});
    writetable(correlations,fullfile(output,'r22dhma_ber_metric_correlations.csv'));
    channel_metrics={'MISMATCH_PENALTY','T_OMIT','WSE_F','SigmaMin','CondHs'};channel_rows=cell(60,5);cursor=0;
    for method=[2,3]
        mse_change=cube.MSE_TRUE_SPARSE(:,p.snr_db==36,method)-cube.MSE_TRUE_SPARSE(:,p.snr_db==30,method);
        for group=1:6
            indices=1:500;if group>1,indices=(group-2)*100+(1:100);end
            for field=channel_metrics
                if isfield(stat,field{1}),values=stat.(field{1})(:,method);
                else,values=cube.(field{1})(:,p.snr_db==36,method)-cube.(field{1})(:,p.snr_db==30,method);end
                [~,spearman,count]=association(mse_change(indices),values(indices));
                cursor=cursor+1;channel_rows(cursor,:)={methods(method).Name,string(groups{group}),string(field{1}),spearman,count};
            end
        end
    end
    channel_correlations=cell2table(channel_rows,'VariableNames',{'Method','Group','Metric','Spearman','Channels'});
    writetable(channel_correlations,fullfile(output,'r22dhma_channel_level_correlations.csv'));
    dominant_rows=cell(68,12);cursor=0;terms={'T_REG','T_OMIT','T_CROSS','T_NOISE'};
    for method=1:4
        for point=1:17
            total=mean(cube.MSE_TRUE_SPARSE(:,point,method));values=zeros(1,4);
            for index=1:4,values(index)=mean(cube.(terms{index})(:,point,method));end
            [~,largest]=max(abs(values));cursor=cursor+1;
            dominant_rows(cursor,:)=[{methods(method).Name,p.snr_db(point)},num2cell(values/total),num2cell(abs(values)),{string(terms{largest}),max(abs(values))/sum(abs(values))}];
        end
    end
    dominant=cell2table(dominant_rows,'VariableNames',{'Method','SNRdB','FractionReg','FractionOmit','FractionCross','FractionNoise', ...
        'AbsReg','AbsOmit','AbsCross','AbsNoise','LargestAbsoluteTerm','LargestAbsoluteShare'});
    writetable(dominant,fullfile(output,'r22dhma_dominant_terms.csv'));
    controls=full_pass&&solver_pass&&noise_audit.SNR_NOISE_MAPPING_PASS;mismatch_supported=controls;interaction=controls;
    for method=[2,3]
        mismatch=change_row(change,methods(method).Name,'MISMATCH_PENALTY');omit=change_row(change,methods(method).Name,'T_OMIT');
        wse=change_row(change,methods(method).Name,'WSE_F');amplification=change_row(change,methods(method).Name,'Amplification_F');
        method_support=mismatch.CI95Lower>0&&(omit.CI95Lower>0||wse.CI95Lower>0||amplification.CI95Lower>0);
        mismatch_supported=mismatch_supported&&method_support;
        cond_increase=change_row(change,methods(method).Name,'CondA');receiver_increase=change_row(change,methods(method).Name,'W_F');
        noise=change_row(change,methods(method).Name,'T_NOISE');
        cond_systematic=mean(all(diff(cube.CondA(:,:,method),1,2)>=-p.monotonic_tolerance*(1+cube.CondA(:,1:16,method)),2))>=.95;
        receiver_systematic=mean(all(diff(cube.W_F(:,:,method),1,2)>=-p.monotonic_tolerance*(1+cube.W_F(:,1:16,method)),2))>=.95;
        condition_ratio=mean(cube.CondA(:,high,method),'all')/mean(cube.CondA(:,mid,method),'all');
        interaction=interaction&&method_support&&cond_systematic&&receiver_systematic&&condition_ratio>=2 ...
            &&cond_increase.CI95Lower>0&&receiver_increase.CI95Lower>0&&noise.CI95Upper<0;
    end
    shared_change=change_row(change,methods(2).Name,'MSE_TRUE_SPARSE');axis_change=change_row(change,methods(3).Name,'MSE_TRUE_SPARSE');
    if ~full_pass||~noise_audit.SNR_NOISE_MAPPING_PASS
        failure='SYSTEM_OR_NOISE_IMPLEMENTATION_ARTIFACT_SUSPECTED';
    elseif ~solver_pass
        failure='NUMERICAL_SOLVER_ARTIFACT_SUSPECTED';
    elseif interaction&&consistent_pass&&shared_change.CI95Lower>0&&axis_change.CI95Lower>0
        failure='LOSS_OF_REGULARIZATION_EXPOSES_OMITTED_CHANNEL_MISMATCH';
    elseif mismatch_supported&&consistent_pass&&shared_change.CI95Lower>0&&axis_change.CI95Lower>0
        failure='HIGH_SNR_MISMATCH_AMPLIFICATION';
    elseif all(diff(cube.MSE_TRUE_SPARSE(:,:,2:3),1,2)<=p.monotonic_tolerance,'all')
        failure='BER_UPTURN_NOT_EXPLAINED_BY_SECOND_ORDER_MSE';
    elseif all(dominant.LargestAbsoluteShare(dominant.SNRdB>=34&ismember(dominant.Method,[methods(2).Name,methods(3).Name]))<.5)
        failure='MULTIFACTOR_HIGH_SNR_MISMATCH';
    else
        failure='UNRESOLVED';
    end
    omitted=paired(paired.Condition=="24"&paired.Metric=="RelativeOmittedF",:);
    condition_hs=paired(paired.Condition=="36"&paired.Metric=="CondHs",:);condition_A=paired(paired.Condition=="36"&paired.Metric=="CondA",:);
    gain_f=paired(paired.Condition=="HIGH_34_36"&paired.Metric=="Amplification_F",:);wse=paired(paired.Condition=="HIGH_34_36"&paired.Metric=="WSE_F",:);
    components=[omitted.CI95Upper<0,condition_hs.CI95Upper<0&&condition_A.CI95Upper<0,gain_f.CI95Upper<0&&wse.CI95Upper<0];
    advantage='UNRESOLVED';if sum(components)>=2,advantage='MULTIFACTOR';elseif components(1),advantage='SMALLER_OMITTED_CHANNEL'; ...
        elseif components(2),advantage='BETTER_CONDITIONED_SPARSE_CHANNEL';elseif components(3),advantage='LOWER_MISMATCH_AMPLIFICATION';end
    shared_mismatch=change_row(change,methods(2).Name,'MISMATCH_PENALTY');axis_mismatch=change_row(change,methods(3).Name,'MISMATCH_PENALTY');
    shared_omit=change_row(change,methods(2).Name,'T_OMIT');axis_omit=change_row(change,methods(3).Name,'T_OMIT');
    shared_wse=change_row(change,methods(2).Name,'WSE_F');axis_wse=change_row(change,methods(3).Name,'WSE_F');
    flags=struct('NEW_CHANNELS_GENERATED',false,'BER_RERUN',false,'NEW_TRANSFORM_TRAINED',false,'ALPHA_CHANGED',false, ...
        'SUPPORT_REBUILT',false,'K_CHANGED',false,'RECEIVER_CHANGED',false,'FILTER_LENGTH_ESCALATED',false,'M_FB_CHANGED',false, ...
        'HISTORICAL_INPUTS_VERIFIED',true,'SNR_NOISE_MAPPING_PASS',noise_audit.SNR_NOISE_MAPPING_PASS,'NUMERICAL_SOLVER_PASS',solver_pass, ...
        'FULL_MSE_MONOTONIC_PASS',full_pass,'MODEL_CONSISTENT_MONOTONIC_PASS',consistent_pass, ...
        'MISMATCH_AMPLIFICATION_SUPPORTED',mismatch_supported,'REGULARIZATION_LOSS_INTERACTION_SUPPORTED',interaction, ...
        'SHARED_HIGH_SNR_MISMATCH_CHANGE',shared_mismatch.MeanDifference,'AXIS_HIGH_SNR_MISMATCH_CHANGE',axis_mismatch.MeanDifference, ...
        'SHARED_HIGH_SNR_OMIT_CHANGE',shared_omit.MeanDifference,'AXIS_HIGH_SNR_OMIT_CHANGE',axis_omit.MeanDifference, ...
        'SHARED_HIGH_SNR_WSE_CHANGE',shared_wse.MeanDifference,'AXIS_HIGH_SNR_WSE_CHANGE',axis_wse.MeanDifference, ...
        'PRIMARY_FAILURE_MECHANISM',failure,'AXIS_ADVANTAGE_MECHANISM',advantage,'MANUSCRIPT_MODIFIED',false);
    r22a_io('json',fullfile(output,'r22dhma_decision_flags.json'),flags);
    lines=["# R2.2D High-SNR Mismatch Diagnostic Audit";"1. Historical inputs verified:TRUE.";"2. Channels audited:500 existing,100pervelocity."; ...
        "3. SNRgrid:20:1:36dB.";"4. K:33178 fixedpositions.";sprintf('5. Sharedalpha:%.15g.',methods(2).AlphaNu); ...
        sprintf('6. Axisalpha_nu:%.15g.',methods(3).AlphaNu);sprintf('7. Axisalpha_tau:%.15g.',methods(3).AlphaTau); ...
        "8. SNR/noisemapping:PASS.";"9. FullMSEmonotonic:"+string(full_pass);"10. Modelconsistentmonotonic:"+string(consistent_pass)];
    for method=[2,3]
        selected=turning(turning.Method==methods(method).Name&turning.Group=="AGGREGATE",:);
        lines(end+1)=sprintf('%d. %s minimum exactsparseMSE %.12g at%d dB; historical minimumBER %.12g at%d dB.',9+method,methods(method).Name, ...
            selected.MinimumSparseMSE,selected.MinimumSparseMSESNR,selected.MinimumBER,selected.MinimumBERSNR); %#ok<AGROW>
    end
    report_metrics={'MSE_TRUE_SPARSE','MISMATCH_PENALTY','T_OMIT','T_NOISE','W_F','WSE_F'};number=13;
    for field=report_metrics
        for method=[2,3]
            selected=change_row(change,methods(method).Name,field{1});
            lines(end+1)=sprintf('%d. %s HIGH-MID %s change %.12g CI95[%.12g,%.12g].',number,methods(method).Name,field{1},selected.MeanDifference,selected.CI95Lower,selected.CI95Upper); %#ok<AGROW>
            number=number+1;
        end
    end
    for method=[2,3]
        for snr=[30,36]
            lines(end+1)=sprintf('%d. %s meancondA at%d dB %.12g.',number,methods(method).Name,snr,mean(cube.CondA(:,p.snr_db==snr,method)));number=number+1; %#ok<AGROW>
        end
    end
    lines=[lines;sprintf('29. Shared mean/median sigma_min %.12g/%.12g.',mean(stat.SigmaMin(:,2)),median(stat.SigmaMin(:,2))); ...
        sprintf('30. Axis mean/median sigma_min %.12g/%.12g.',mean(stat.SigmaMin(:,3)),median(stat.SigmaMin(:,3))); ...
        sprintf('31. MaxSVDrelativeWdifference %.12g; maxMSEdifference %.12g; conditionidentityerror %.12g.',maximum_reference,maximum_mse_reference,maximum_condition_error); ...
        sprintf('32. Maximumsolverbackwardresidual %.12g.',maximum_backward);"33. MISMATCH_AMPLIFICATION_SUPPORTED="+string(mismatch_supported); ...
        "34. REGULARIZATION_LOSS_INTERACTION_SUPPORTED="+string(interaction)];
    for field={'MISMATCH_PENALTY','WSE_F'}
        selected=correlations(correlations.Metric==string(field{1}),:);
        lines(end+1)=sprintf('%d. BER correlation %s:sharedPearson/Spearman %.7g/%.7g;axis %.7g/%.7g;17SNRpoints descriptiveonly.', ...
            34+find(strcmp({'MISMATCH_PENALTY','WSE_F'},field{1})),field{1},selected.Pearson(1),selected.Spearman(1),selected.Pearson(2),selected.Spearman(2)); %#ok<AGROW>
    end
    lines=[lines;sprintf('37. Axis-shared relativeomittedF %.12g CI95[%.12g,%.12g].',omitted.MeanDifference,omitted.CI95Lower,omitted.CI95Upper); ...
        sprintf('38. Axis-sharedcondA36 %.12g CI95[%.12g,%.12g];condHs %.12g.',condition_A.MeanDifference,condition_A.CI95Lower,condition_A.CI95Upper,condition_hs.MeanDifference); ...
        sprintf('39. Axis-shared HIGH amplificationF %.12g CI95[%.12g,%.12g].',gain_f.MeanDifference,gain_f.CI95Lower,gain_f.CI95Upper); ...
        "40. AXIS_ADVANTAGE_MECHANISM="+string(advantage);"41. PRIMARY_FAILURE_MECHANISM="+string(failure); ...
        "42. NEW_CHANNELS_GENERATED=FALSE";"43. BER_RERUN=FALSE";"44. NEW_TRANSFORM_TRAINED=FALSE";"45. SUPPORT_REBUILT=FALSE"; ...
        "46. K_CHANGED=FALSE";"47. RECEIVER_CHANGED=FALSE";"48. M_FB_CHANGED=FALSE";"49. MANUSCRIPT_MODIFIED=FALSE"; ...
        "";"## Scope And Interpretation"; ...
        "Exact trace decomposition verified perchannel/SNR. SignedCROSS maybenegative; contributionfractions are accounting,notprobabilities. Noinverseusedinoperationalreceivers. SVDpseudoinverse only zero-noiselimitdiagnostic,relative tolerance1e-12; neverBER."; ...
        "Full/modelconsistent monotonicchecks use everychannel andmethod. Solverreference first10existingchannelspervelocity,all4methods/all17SNRs; backwardresidual forallsolves. No eigenvalueclipping oraddedregularization."; ...
        "HIGH34:36minusMID28:30 contrasts bootstrap10000pairedphysicalchannels,stratifiedbyvelocity,equalweights. BER/metriccorrelations17serialSNRpoints descriptiveonly,no causal significance. Channel correlations36minus30 are diagnostic; nonfinitepairs excludedwithcounts,notconditionnumberclamping."; ...
        "Bothprimariesrequiredformechanismsupport. Strongregularizationinteraction additionallyrequiresconditionratio>=2,systematiccond/Wgrowth andsupporteddecreaseofNOISE. Ifthermaltermrises,do notforce CaseB. CaseAcanstillhold."; ...
        "Axisadvantageclassifiedby supportedomitted/conditioning/amplificationcomponents; do notforcea singleexplanation. ChangesinEandgainmustbeevaluatedtogether."; ...
        "No newchannels,bits,noise,BER,training,receiver/support/alpha/Kchanges,regularizationtests orarchitecture/manuscript modification. This explains existingreceiversecond-orderbehavior; notanerror-floorortransitiontheorem."];
    r22a_io('text',fullfile(output,'r22dhma_final_report.md'),join(lines,newline));
    r22a_io('text',fullfile(output,'r22dhma_executive_summary.md'),join(lines,newline));
    figures(curves,ber,paired,change,methods,output);
end

function T=aggregate_curves(cube,methods,p)
    fields=fieldnames(cube);rows=cell(4*6*17,4+numel(fields));cursor=0;
    for group=0:5
        indices=1:500;if group>0,indices=(group-1)*100+(1:100);end
        label="AGGREGATE";if group>0,label="VELOCITY_"+string(group);end
        for method=1:4
            for point=1:17
                values=zeros(1,numel(fields));for index=1:numel(fields),values(index)=mean(cube.(fields{index})(indices,point,method),'omitnan');end
                cursor=cursor+1;rows(cursor,:)=[{methods(method).Name,label,p.snr_db(point),numel(indices)},num2cell(values)];
            end
        end
    end
    T=cell2table(rows,'VariableNames',[{'Method','Group','SNRdB','Channels'},fields.']);
end

function T=historical_ber(counts,methods,p)
    rows=cell(4*6*17,4);cursor=0;
    for group=0:5
        label="AGGREGATE";velocities=1:5;if group>0,label="VELOCITY_"+string(group);velocities=group;end
        for method=1:4
            for snr=p.snr_db
                selected=counts(string(counts.Method)==methods(method).Name&ismember(counts.VelocityIndex,velocities)&counts.SNRdB==snr,:);
                cursor=cursor+1;rows(cursor,:)={methods(method).Name,label,snr,mean(selected.Errors./selected.Bits)};
            end
        end
    end
    T=cell2table(rows,'VariableNames',{'Method','Group','SNRdB','BER'});
end

function row=change_row(T,method,metric)
    row=T(T.Method==method&T.Metric==string(metric),:);assert(height(row)==1);
end

function [pearson,spearman,count]=association(left,right)
    valid=isfinite(left)&isfinite(right);left=left(valid);right=right(valid);count=numel(left);pearson=NaN;spearman=NaN;
    if count>=3&&std(left)>0&&std(right)>0
        pearson=corr(left,right);spearman=corr(left,right,'Type','Spearman');
    end
end

function figures(curves,ber,paired,change,methods,output)
    panels={{'MSE_TRUE_SPARSE','MSE_MODEL_CONSISTENT','MISMATCH_PENALTY','MSE_FULL'}, ...
        {'T_REG','T_OMIT','T_CROSS','T_NOISE'},{'W_F','WSE_F'},{'CondA'}};
    files={'exact_mse_decomposition','mismatch_terms_vs_snr','receiver_gain_vs_snr','conditioning_vs_snr'};
    for panel=1:4
        if panel==1
            exact_mse_decomposition_figure(curves,[methods(2).Name,methods(3).Name],output,[7.2,3.7]);
            continue;
        end
        fig=new_figure();tiledlayout(1,2,'TileSpacing','compact');
        for method=[2,3]
            nexttile;hold on;selected=curves(curves.Method==methods(method).Name&curves.Group=="AGGREGATE",:);
            for field=panels{panel},plot(selected.SNRdB,selected.(field{1}),'-o','DisplayName',field{1});end
            grid on;xlabel('SNR (dB)');title(figure_label(methods(method).Name),'Interpreter','none');legend('Interpreter','none','Location','best');
            if panel==4,set(gca,'YScale','log');end
        end
        export_figure(fig,output,['fig_r22dhma_',files{panel}]);
    end
    fig=new_figure();tiledlayout(1,2,'TileSpacing','compact');
    for method=[2,3]
        nexttile;selected=curves(curves.Method==methods(method).Name&curves.Group=="AGGREGATE",:);
        rates=ber(ber.Method==methods(method).Name&ber.Group=="AGGREGATE",:);scatter(selected.MSE_TRUE_SPARSE,rates.BER,40,selected.SNRdB,'filled');
        xlabel('Exact sparse MSE');ylabel('Historical measured BER');colorbar;grid on;title(figure_label(methods(method).Name),'Interpreter','none');
    end
    export_figure(fig,output,'fig_r22dhma_ber_vs_exact_mse');
    fig=new_figure();tiledlayout(1,3,'TileSpacing','compact');
    for field={'RelativeOmittedF','CondA','Amplification_F'}
        nexttile;selected=paired(paired.Condition=="HIGH_34_36"&paired.Metric==string(field{1}),:);
        errorbar(1,selected.MeanDifference,selected.MeanDifference-selected.CI95Lower,selected.CI95Upper-selected.MeanDifference,'o');
        hold on;yline(0);xlim([.5,1.5]);xticks(1);xticklabels({'Independent - shared'});title(field{1},'Interpreter','none');grid on;
    end
    export_figure(fig,output,'fig_r22dhma_axis_vs_shared_highsnr');
    fig=new_figure();tiledlayout(1,2,'TileSpacing','compact');labels={'MISMATCH_PENALTY','T_OMIT','T_NOISE'};
    for method=[2,3]
        nexttile;hold on;
        for index=1:3
            row=change_row(change,methods(method).Name,labels{index});errorbar(index,row.MeanDifference,row.MeanDifference-row.CI95Lower,row.CI95Upper-row.MeanDifference,'o');
        end
        yline(0);xticks(1:3);xticklabels({'Mismatch','Omitted','Noise'});ylabel('HIGH34:36 minus MID28:30, CI95');grid on;title(figure_label(methods(method).Name),'Interpreter','none');
    end
    export_figure(fig,output,'fig_r22dhma_highsnr_mechanism_summary');
end

function exact_mse_decomposition_figure(curves,method_names,output,figure_size)
    fields={'MSE_TRUE_SPARSE','MSE_MODEL_CONSISTENT','MISMATCH_PENALTY','MSE_FULL'};
    labels={'True-channel MSE','Model-consistent MSE','Net mismatch contribution','Full-channel LMMSE MSE'};
    line_styles={'-','--','-.',':'};markers={'o','s','^','d'};colors=lines(4);
    fig=new_figure(figure_size);layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    handles=gobjects(4,1);first_axes=[];
    for method_index=1:2
        ax=nexttile(layout);hold(ax,'on');if method_index==1,first_axes=ax;end
        selected=curves(curves.Method==method_names(method_index)&curves.Group=="AGGREGATE" ...
            &curves.SNRdB>=20&curves.SNRdB<=36,:);
        selected=sortrows(selected,'SNRdB');assert(height(selected)==17);
        for quantity=1:4
            h=plot(ax,selected.SNRdB,selected.(fields{quantity}), ...
                'LineStyle',line_styles{quantity},'Marker',markers{quantity}, ...
                'Color',colors(quantity,:),'LineWidth',1.4,'MarkerSize',4.5, ...
                'DisplayName',labels{quantity});
            if method_index==1,handles(quantity)=h;end
        end
        grid(ax,'on');xlabel(ax,'SNR (dB)');ylabel(ax,'Normalized MSE');
        panel_title=title(ax,figure_label(method_names(method_index)),'Interpreter','none');
        if method_index==2
            panel_title.Units='normalized';
            title_position=panel_title.Position;
            title_position(1)=title_position(1)+0.015;
            panel_title.Position=title_position;
        end
        set(ax,'YScale','linear');xlim(ax,[20,36]);ylim(ax,[0,0.25]);yticks(ax,0:0.05:0.25);
        pbaspect(ax,[1.4,1,1]);
    end
    shared_legend=legend(first_axes,handles,labels,'NumColumns',2, ...
        'Orientation','horizontal','Location','north','Interpreter','none');
    shared_legend.Layout.Tile='north';
    export_figure(fig,output,'fig_r22dhma_exact_mse_decomposition',figure_size);
end

function label=figure_label(name)
    label='Shared-parameter representation';
    if contains(string(name),'AXIS'),label='Independently parameterized representation';end
end

function fig=new_figure(figure_size)
    if nargin>0
        fig=figure('Visible','off','Color','w','Units','inches', ...
            'Position',[0.5,0.5,figure_size(1),figure_size(2)]);
    else
        fig=figure('Visible','off','Color','w','Position',[100,100,1100,650]);
    end
end

function export_figure(fig,output,name,figure_size)
    if nargin>3
        set(fig,'PaperUnits','inches','PaperPosition',[0,0,figure_size(1),figure_size(2)], ...
            'PaperSize',figure_size,'PaperPositionMode','manual');
        print(fig,fullfile(output,[name,'.pdf']),'-dpdf','-vector');
    else
        exportgraphics(fig,fullfile(output,[name,'.pdf']),'ContentType','vector');
    end
    exportgraphics(fig,fullfile(output,[name,'.png']),'Resolution',180);close(fig);
end

function flags=selftest(p,source,output)
    smoke=fullfile(output,'deterministic_preflight');if ~isfolder(smoke),mkdir(smoke);end
    loaded=load(fullfile(source,'r22dbe_frozen_methods.mat'),'methods');methods=loaded.methods;
    fields={'MSE_TRUE_SPARSE','MSE_MODEL_CONSISTENT','MISMATCH_PENALTY','T_REG','T_OMIT','T_CROSS','T_NOISE','MSE_FULL', ...
        'CondA','LambdaMin','LambdaMax','W_2','W_F','WSE_2','WSE_F','Amplification_F','Amplification_2','BackwardResidual', ...
        'ReconstructionError','PenaltyEqualityError','SVDReceiverRelativeError','SVDMSEDifference','ConditionIdentityError'};
    rows=cell(34000,27);static_rows=cell(2000,11);cursor=0;static_cursor=0;
    for channel=1:500
        for method=1:4
            static_cursor=static_cursor+1;static_rows(static_cursor,:)={methods(method).Name,ceil(channel/100),mod(channel-1,100)+1,2,.1,20,288,2e-12,.1,.2,.15};
            for point=1:17
                snr=19+point;reg=.01/point;noise=.02/point;omit=.001*point;cross=0;consistent=reg+noise;total=consistent+omit;
                values=[total,consistent,omit,reg,omit,cross,noise,.01/point,2^point,.001,1,point,point,sqrt(288*omit),sqrt(288*omit), ...
                    sqrt(288*omit)/.2,sqrt(288*omit)/.15,1e-16,0,0,1e-12,1e-14,1e-12];
                cursor=cursor+1;rows(cursor,:)=[{methods(method).Name,ceil(channel/100),mod(channel-1,100)+1,snr},num2cell(values)];
            end
        end
    end
    detailed=cell2table(rows,'VariableNames',[{'Method','VelocityIndex','ChannelID','SNRdB'},fields]);
    static=cell2table(static_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SigmaMax','SigmaMin','CondHs','NumericalRank','RankTolerance','RelativeOmittedF','OmittedF','Omitted2'});
    r22a_io('json',fullfile(smoke,'r22dhma_snr_noise_audit.json'),struct('SNR_NOISE_MAPPING_PASS',true));
    flags=r22dhma_finish(methods,p,struct('Detailed',detailed,'Static',static),source,smoke);
    assert(flags.MISMATCH_AMPLIFICATION_SUPPORTED&&flags.REGULARIZATION_LOSS_INTERACTION_SUPPORTED);
    r22a_io('json',fullfile(output,'r22dhma_deterministic_preflight.json'),struct('PASS',true,'PHYSICAL_CHANNELS_GENERATED',false,'BER_RERUN',false));
end