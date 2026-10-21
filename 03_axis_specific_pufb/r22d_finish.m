function flags=r22d_finish(methods,p,mse,frob,sparsity,full,stability,output)
    if ischar(methods)&&strcmp(methods,'SELFTEST')
        flags=synthetic_check(p,output);return
    end
    delta=zeros(1000,8,6);full_values=delta;frob_values=zeros(1000,6);k99=frob_values;j4=frob_values;coverage=frob_values;
    for method=1:6
        name=methods(method).Name;selected=sortrows(mse(string(mse.Method)==name,:),{'VelocityIndex','ChannelID','SNRdB'});
        assert(height(selected)==8000&&all(selected.RetainedEntries==p.K));
        assert(isequal(selected.VelocityIndex,repelem((1:5)',1600))&&isequal(selected.ChannelID,repmat(repelem((1:200)',8),5,1)) ...
            &&isequal(selected.SNRdB,repmat((3:3:24)',1000,1)),'R22D:Pairing','Incomplete or duplicated TEST keys.');
        delta(:,:,method)=reshape(selected.DeltaMSE,8,1000).';
        selected=sortrows(full(string(full.Method)==name,:),{'VelocityIndex','ChannelID','SNRdB'});
        full_values(:,:,method)=reshape(selected.FullMSE,8,1000).';
        selected_frob=sortrows(frob(string(frob.Method)==name,:),{'VelocityIndex','ChannelID'});
        selected_sparse=sortrows(sparsity(string(sparsity.Method)==name,:),{'VelocityIndex','ChannelID'});
        assert(height(selected_frob)==1000&&height(selected_sparse)==1000);
        frob_values(:,method)=selected_frob.FrobeniusError;k99(:,method)=selected_sparse.K99;
        j4(:,method)=selected_sparse.J4;coverage(:,method)=selected_sparse.SupportCoverage;
    end
    invariance=max(abs(full_values-full_values(:,:,1)),[],'all');assert(invariance<p.full_invariance_tolerance);
    groups=struct('Label',{'AGGREGATE','LOW','MID','HIGH'},'Velocities',{1:5,1:5,1:5,1:5},'Points',{[],1:3,4:6,7:8});
    for velocity=1:5,groups(end+1)=struct('Label',sprintf('VELOCITY_%d',velocity),'Velocities',velocity,'Points',[]);end %#ok<AGROW>
    groups(end+1)=struct('Label','MACH5_21_24','Velocities',5,'Points',7:8);
    for point=1:8,groups(end+1)=struct('Label',sprintf('SNR_%d',point*3),'Velocities',1:5,'Points',point);end %#ok<AGROW>
    rows=cell(5*numel(groups),12);cursor=0;
    for reference=[4,6,1,2,3]
        for group=1:numel(groups)
            points=groups(group).Points;velocities=groups(group).Velocities;
            left=channel_average(delta(:,:,reference),points);right=channel_average(delta(:,:,5),points);
            indices=reshape((1:1000)',200,5);indices=indices(:,velocities);indices=indices(:);
            difference=right(indices)-left(indices);[lower,upper]=bootstrap(difference,numel(velocities),p);
            cursor=cursor+1;rows(cursor,:)={"AXIS_MINUS_"+methods(reference).Name,string(groups(group).Label), ...
                mean(difference),lower,upper,mean(left(indices)),mean(right(indices)), ...
                (mean(left(indices))-mean(right(indices)))/mean(left(indices)),median(difference),std(difference), ...
                mean(difference)/max(std(difference),eps),numel(difference)};
        end
    end
    paired=cell2table(rows,'VariableNames',{'Comparison','Group','MeanDifference','CI95Lower','CI95Upper','ReferenceMean','CandidateMean', ...
        'RelativeReduction','PairedMedian','PairedSD','StandardizedEffect','PhysicalChannels'});
    assert(all(isfinite(paired.MeanDifference))&&all(isfinite(paired.CI95Lower))&&all(isfinite(paired.CI95Upper)));
    writetable(paired,fullfile(output,'r22d_pairwise_bootstrap.csv'));
    writetable(paired(startsWith(paired.Group,"SNR_"),:),fullfile(output,'r22d_snr_breakdown.csv'));
    writetable(paired(startsWith(paired.Group,"VELOCITY_"),:),fullfile(output,'r22d_velocity_breakdown.csv'));
    writetable(paired(paired.Group=="MACH5_21_24",:),fullfile(output,'r22d_stress_cell.csv'));
    primary=contrast(paired,methods(4).Name,'AGGREGATE');swap=contrast(paired,methods(6).Name,'AGGREGATE');
    stress=contrast(paired,methods(4).Name,'MACH5_21_24');historical=contrast(paired,methods(3).Name,'AGGREGATE');
    separation=abs(methods(5).AlphaTau-methods(5).AlphaNu);circular=min(separation,pi-separation);
    separated=separation>p.parameter_separation_multiple*p.final_axis_mesh;
    benefit=classify(primary);assignment=classify(swap);stress_benefit=classify(stress);
    specialization=strcmp(benefit,'SUPPORTED')&&separated&&strcmp(assignment,'SUPPORTED');
    secondary=paired(paired.Comparison=="AXIS_MINUS_"+methods(4).Name& ...
        (ismember(paired.Group,["LOW","MID","HIGH"])|startsWith(paired.Group,"VELOCITY_")),:);
    regime=strcmp(benefit,'SUPPORTED')&&any(secondary.CI95Lower>=0);
    if regime
        interpretation='AXIS_SPECIFIC_RECEIVER_AWARE_GAIN_IS_REGIME_DEPENDENT';
    elseif specialization
        interpretation='DELAY_AND_DOPPLER_BENEFIT_FROM_DISTINCT_RECEIVER_AWARE_PUFB_GEOMETRIES';
    elseif strcmp(benefit,'SUPPORTED')&&~separated
        interpretation='AXIS_SPECIFIC_PARAMETERIZATION_IMPROVES_NUMERIC_SEARCH_WITHOUT_RESOLVED_PARAMETER_SEPARATION';
    elseif strcmp(benefit,'SUPPORTED')
        interpretation='SECOND_DEGREE_OF_FREEDOM_IMPROVES_RECEIVER_ERROR_WITHOUT_CONFIRMED_AXIS_SPECIFICITY';
    elseif circular<=p.final_axis_mesh
        interpretation='RECEIVER_OPTIMUM_REMAINS_EFFECTIVELY_SHARED_ACROSS_AXES';
    elseif strcmp(classify(historical),'SUPPORTED')
        interpretation='NO_AXIS_SPECIFIC_GAIN_AFTER_CONTROLLING_FOR_DATA_AND_REOPTIMIZATION';
    else
        interpretation='SHARED_MOTHER_BANK_IS_SUFFICIENT_WITHIN_THE_LW4_RECEIVER_AWARE_FAMILY';
    end
    next='NONE_PENDING_REASSESSMENT';if strcmp(benefit,'SUPPORTED'),next='R2.2D_CONFIRM';end
    flags=struct('PHYSICAL_SCENARIO','ABID_EVA','FILTER_LENGTH',4,'J_NU',2,'J_TAU',3,'PRIMARY_RETAINED_ENTRIES',p.K, ...
        'SHARED_CONTROL_TRAINED_ON_SAME_BANK',true,'AXIS_SPECIFIC_PARAMETERS',true, ...
        'ALPHA_NU_STAR',methods(5).AlphaNu,'ALPHA_TAU_STAR',methods(5).AlphaTau,'ALPHA_SHARED_STAR',methods(4).AlphaNu, ...
        'ABS_ALPHA_SEPARATION',separation,'CIRCULAR_ALPHA_SEPARATION',circular,'PARAMETERS_NUMERICALLY_SEPARATED',separated, ...
        'AXIS_SPECIFIC_MSE_BENEFIT',benefit,'AXIS_ASSIGNMENT_MATTERS',assignment,'AXIS_SPECIALIZATION_SUPPORTED',specialization, ...
        'AXIS_SPECIFIC_STRESS_BENEFIT',stress_benefit,'FULL_RECEIVER_INVARIANCE_PASS',invariance<1e-12, ...
        'FULL_RECEIVER_INVARIANCE_RESIDUAL',invariance,'BER_RUN',false,'PCG_RUN',false,'K99_USED_FOR_SELECTION',false, ...
        'BER_USED_FOR_TRAINING',false,'GRAM_PENALTY_USED',false,'FILTER_LENGTH_ESCALATED',false, ...
        'NONSEPARABLE_TRANSFORM_USED',false,'TEST_GENERATED_AFTER_COMPLETE_FREEZE',true,'TEST_USED_FOR_SELECTION',false, ...
        'NEXT_EXPERIMENT',next,'MANUSCRIPT_MODIFIED',false,'PRIMARY_INTERPRETATION',interpretation);
    r22a_io('json',fullfile(output,'r22d_decision_flags.json'),flags);
    training_shared=r22d_model('read',fullfile(output,'r22d_shared_coarse_search.csv'));
    training_axis=r22d_model('read',fullfile(output,'r22d_axis_global_surface.csv'));
    refined=r22d_model('read',fullfile(output,'r22d_axis_refinement.csv'));
    selection=r22d_model('read',fullfile(output,'r22d_val_selection.csv'));
    shared_selection=sortrows(selection(selection.Arm=="SHARED",:),{'ValJrx','AlphaNu','AlphaTau'});
    axis_selection=sortrows(selection(selection.Arm=="AXIS",:),{'ValJrx','AlphaNu','AlphaTau'});
    means=zeros(1,6);for method=1:6,means(method)=mean(channel_average(delta(:,:,method),[]));end
    lines=["# R2.2D Axis-Specific Receiver-Aware Fixed-Budget Pilot"; ...
        "1. Canonical alpha interval independently on each axis: [0,pi)."; ...
        "2. TRAIN_SCREEN size: 100 (20 per velocity).";"3. Full TRAIN size: 200 (40 per velocity)."; ...
        "4. VAL size: 100 (20 per velocity), one selection access for both arms.";"5. TEST size: 1000 (200 per velocity)."; ...
        "6. K: exactly 33178 common TRAIN coefficients per method."; ...
        sprintf('7. Shared coarse candidates: %d.',height(training_shared)); ...
        sprintf('8. Axis-specific global candidates including mandatory insertions: %d.',height(training_axis)); ...
        sprintf('9. Distinct axis-specific refined candidates: %d.',height(refined)); ...
        sprintf('10. Same-bank shared alpha: %.15g.',methods(4).AlphaNu); ...
        sprintf('11. Shared VAL J_RX: %.12g.',shared_selection.ValJrx(1)); ...
        sprintf('12. alpha_nu_star: %.15g.',methods(5).AlphaNu); ...
        sprintf('13. alpha_tau_star: %.15g.',methods(5).AlphaTau); ...
        sprintf('14. |alpha_nu-alpha_tau|: %.12g; signed tau-minus-nu: %.12g; final mesh %.12g, separation threshold %.12g.', ...
            separation,methods(5).AlphaTau-methods(5).AlphaNu,p.final_axis_mesh,p.parameter_separation_multiple*p.final_axis_mesh); ...
        sprintf('15. Axis-specific VAL J_RX: %.12g.',axis_selection.ValJrx(1)); ...
        sprintf('16. Shared lowpass filter: [%.12g %.12g %.12g %.12g].',methods(4).Nu.h); ...
        sprintf('17. Doppler lowpass filter: [%.12g %.12g %.12g %.12g].',methods(5).Nu.h); ...
        sprintf('18. Delay lowpass filter: [%.12g %.12g %.12g %.12g].',methods(5).Tau.h); ...
        sprintf('19. Axis paraunitarity nu/tau: %.12g / %.12g.',methods(5).Nu.ParaunitarityResidual,methods(5).Tau.ParaunitarityResidual); ...
        sprintf('20. Axis PR residual ||UU''-I||_F: %.12g; max structural gate %.12g.',methods(5).Residuals(4),max(methods(5).Residuals)); ...
        sprintf('21. D4 TEST DeltaMSE: %.12g.',means(1));sprintf('22. Historical R2.2A on NEW TEST: %.12g.',means(2)); ...
        sprintf('23. Confirmed R2.2C on NEW TEST: %.12g.',means(3));sprintf('24. R22D shared TEST: %.12g.',means(4)); ...
        sprintf('25. R22D axis-specific TEST: %.12g.',means(5));sprintf('26. Axis-minus-shared difference: %.12g.',primary.MeanDifference); ...
        sprintf('27. Relative reduction: %.9g%%.',100*primary.RelativeReduction); ...
        sprintf('28. Paired CI95: [%.12g, %.12g].',primary.CI95Lower,primary.CI95Upper); ...
        sprintf('29. Standardized paired effect: %.12g; median %.12g; sample SD %.12g.',primary.StandardizedEffect,primary.PairedMedian,primary.PairedSD)];
    for index=1:3
        labels={'LOW','MID','HIGH'};selected=contrast(paired,methods(4).Name,labels{index});
        lines(end+1)=sprintf('%d. %s SNR difference %.12g CI95 [%.12g,%.12g].',29+index,labels{index},selected.MeanDifference,selected.CI95Lower,selected.CI95Upper); %#ok<AGROW>
    end
    for velocity=1:5
        selected=contrast(paired,methods(4).Name,sprintf('VELOCITY_%d',velocity));
        lines(end+1)=sprintf('%d. %s difference %.12g CI95 [%.12g,%.12g].',32+velocity,p.velocity_labels{velocity},selected.MeanDifference,selected.CI95Lower,selected.CI95Upper); %#ok<AGROW>
    end
    shared_stability=stability(string(stability.Method)==methods(4).Name,:);axis_stability=stability(string(stability.Method)==methods(5).Name,:);
    lines=[lines;sprintf('38. Mach5 21/24 difference %.12g CI95 [%.12g,%.12g].',stress.MeanDifference,stress.CI95Lower,stress.CI95Upper); ...
        sprintf('39. Swapped-axis TEST DeltaMSE: %.12g.',means(6)); ...
        sprintf('40. STAR-minus-SWAP %.12g CI95 [%.12g,%.12g].',swap.MeanDifference,swap.CI95Lower,swap.CI95Upper); ...
        sprintf('41. K99 shared/axis: %.12g / %.12g.',mean(k99(:,4)),mean(k99(:,5))); ...
        sprintf('42. Frobenius error shared/axis: %.12g / %.12g.',mean(frob_values(:,4)),mean(frob_values(:,5))); ...
        sprintf('43. J4 shared/axis: %.12g / %.12g.',mean(j4(:,4)),mean(j4(:,5))); ...
        sprintf('44. Support coverage shared/axis %.12g/%.12g; within Jaccard %.12g/%.12g; cross Jaccard %.12g/%.12g.', ...
            mean(coverage(:,4)),mean(coverage(:,5)),mean(shared_stability.WithinVelocityJaccard),mean(axis_stability.WithinVelocityJaccard), ...
            mean(shared_stability.CrossVelocityJaccard),mean(axis_stability.CrossVelocityJaccard)); ...
        sprintf('45. Full-receiver invariance residual %.12g.',invariance); ...
        "46. PARAMETERS_NUMERICALLY_SEPARATED="+string(separated);"47. AXIS_SPECIFIC_MSE_BENEFIT="+string(benefit); ...
        "48. AXIS_ASSIGNMENT_MATTERS="+string(assignment);"49. AXIS_SPECIALIZATION_SUPPORTED="+string(specialization); ...
        "50. AXIS_SPECIFIC_STRESS_BENEFIT="+string(stress_benefit);"51. PRIMARY_INTERPRETATION="+string(interpretation); ...
        "52. NEXT_EXPERIMENT="+string(next);"53. BER_RUN=FALSE";"54. PCG_RUN=FALSE";"55. FILTER_LENGTH_ESCALATED=FALSE"; ...
        "56. TEST_USED_FOR_SELECTION=FALSE";"57. MANUSCRIPT_MODIFIED=FALSE"; ...
        "";"## Statistical And Procedural Notes"; ...
        "Primary mean equally weights LOW/MID/HIGH regimes within each physical channel and five fixed velocities. Physical channel is the paired bootstrap unit; 10000 stratified resamples. SNR points from one channel are not independent units. Secondary intervals are unadjusted/descriptive."; ...
        "Only shared vs axis on the same new TRAIN/VAL/TEST banks is the primary comparison. All masks use the same K and full-TRAIN construction; historical alphas are immutable secondary comparators. SWAP was built after both arms were frozen, not used in selection."; ...
        "Numerical separation alone does not establish physical specialization: benefit, separation and assignment test must all pass. Null results imply sufficiency only within this family/search and power, not proof of identical physical geometries."; ...
        "Surface heatmap/contours show the100-channel/3SNR coarse screen. Full200-channel/8SNR refinement samples are separately marked, not conflated with screen scores. Rendering interpolation is visualization only."; ...
        "Jaccards use first40 TEST channels per velocity, preregistered for tractable descriptive diagnostics. Oracle supports are diagnostic only and never used by receivers or selection."; ...
        sprintf('Coarse axis spacing %.12g; shared final mesh %.12g; axis final mesh %.12g.',pi/17,pi/(61*10),p.final_axis_mesh); ...
        "Historical shared-alpha visual0.1%regret interval is [0.783701703365,0.789676909589], width0.005975206224. Compare descriptive scales only, not a separation hypothesis threshold."; ...
        "No BER/PCG, manuscript edits, filter escalation, black-box optimizer or additional parameters. Next experiment is recommendation only."];
    r22a_io('text',fullfile(output,'r22d_final_report.md'),join(lines,newline));
    r22a_io('text',fullfile(output,'r22d_executive_summary.md'),join(lines,newline));
    figures(methods,p,delta,paired,frob_values,k99,output,selection);
end

function values=channel_average(matrix,points)
    if isempty(points),values=(mean(matrix(:,1:3),2)+mean(matrix(:,4:6),2)+mean(matrix(:,7:8),2))/3;
    else,values=mean(matrix(:,points),2);end
end

function [lower,upper]=bootstrap(values,velocity_count,p)
    values=reshape(values,200,velocity_count);stream=RandStream('mt19937ar','Seed',p.seeds.BOOTSTRAP);
    means=zeros(p.bootstrap_replicates,velocity_count);
    for velocity=1:velocity_count
        draws=randi(stream,200,[200,p.bootstrap_replicates]);current=values(:,velocity);means(:,velocity)=mean(current(draws),1).';
    end
    samples=sort(mean(means,2));lower=samples(ceil(.025*p.bootstrap_replicates));upper=samples(ceil(.975*p.bootstrap_replicates));
end

function row=contrast(T,reference,group)
    row=T(T.Comparison=="AXIS_MINUS_"+string(reference)&T.Group==string(group),:);assert(height(row)==1);
end

function value=classify(row)
    value='INCONCLUSIVE';if row.CI95Upper<0,value='SUPPORTED';elseif row.CI95Lower>=0,value='NOT_SUPPORTED';end
end

function figures(methods,p,delta,paired,frob_values,k99,output,selection)
    surface=r22d_model('read',fullfile(output,'r22d_axis_global_surface.csv'));coarse_grid=surface(1:289,:);
    axis_values=(0:16)*pi/17;values=reshape(coarse_grid.Jrx,17,17).';
    refined=r22d_model('read',fullfile(output,'r22d_axis_refinement.csv'));[~,best]=min(refined.Jrx);
    shared_train=r22d_model('read',fullfile(output,'r22d_shared_refinement.csv'));[~,best_shared]=min(shared_train.Jrx);
    points=[p.historical_r22c_alpha,p.historical_r22c_alpha;shared_train.AlphaNu(best_shared),shared_train.AlphaTau(best_shared);refined.AlphaNu(best),refined.AlphaTau(best)];
    labels={'Historical shared','Best shared TRAIN','Best independent TRAIN'};
    fig=new_figure();imagesc(axis_values,axis_values,values);axis xy;axis square;colorbar;hold on;
    plot([0,pi],[0,pi],'w--','DisplayName','Shared diagonal');overlay(points,labels);xlabel('alpha_nu');ylabel('alpha_tau');title('TRAIN_SCREEN J_{RX}, 100 channels / 3 SNRs');legend('Location','best');
    export_figure(fig,output,'fig_r22d_train_surface_heatmap');
    fig=new_figure();contour(axis_values,axis_values,values,20,'ShowText','on');hold on;plot([0,pi],[0,pi],'k--','DisplayName','Shared diagonal');
    overlay(points,labels);xlabel('alpha_nu');ylabel('alpha_tau');title('TRAIN_SCREEN contour map');grid on;legend('Location','best');
    export_figure(fig,output,'fig_r22d_train_surface_contours');
    fig=new_figure();surf(axis_values,axis_values,values,'EdgeColor','none','FaceAlpha',.75,'DisplayName','Coarse TRAIN_SCREEN');hold on;
    scatter3(refined.AlphaNu,refined.AlphaTau,refined.Jrx,12,'filled','DisplayName','Full-TRAIN refinement');
    diagonal_values=interp2(axis_values,axis_values,values,axis_values,axis_values);
    plot3(axis_values,axis_values,diagonal_values,'k--','DisplayName','Shared screen diagonal');
    for point=1:3
        plot_height=interp2(axis_values,axis_values,values,points(point,1),points(point,2));
        scatter3(points(point,1),points(point,2),plot_height,70,'filled','DisplayName',labels{point});
    end
    xlabel('alpha_nu');ylabel('alpha_tau');zlabel('Stage-specific J_{RX}');title('Screen surface and separately marked full-TRAIN samples');view(40,30);grid on;legend('Location','best');
    export_figure(fig,output,'fig_r22d_train_surface_3d');
    fig=new_figure();bar(1:height(selection),selection.ValJrx);xticks(1:height(selection));
    xticklabels(compose('%s %d',replace(selection.Arm,["AXIS","SHARED"],["Independent","Shared"]),selection.CandidateID));ylabel('VAL J_{RX}');grid on;
    export_figure(fig,output,'fig_r22d_shared_vs_axis_val');
    colors=[.1,.3,.6;.65,.22,.13;.2,.55,.3];primary_methods=[4,5,6];labels={'Shared','Independent','Parameter-swapped'};
    fig=new_figure();tiledlayout(2,3,'TileSpacing','compact');
    for velocity=1:5
        nexttile;hold on;indices=(velocity-1)*200+(1:200);
        for index=1:3,semilogy(p.snr_db,mean(delta(indices,:,primary_methods(index)),1),'-o','Color',colors(index,:),'DisplayName',labels{index});end
        set(gca,'YScale','log');xlabel('SNR (dB)');ylabel('Exact excess MSE');title(p.velocity_labels{velocity});grid on;
        if velocity==1,legend('Location','best');end
    end
    export_figure(fig,output,'fig_r22d_test_excess_mse_vs_snr');
    fig=new_figure();selected=paired(paired.Comparison=="AXIS_MINUS_"+methods(4).Name&startsWith(paired.Group,"SNR_"),:);
    errorbar(p.snr_db,selected.MeanDifference,selected.MeanDifference-selected.CI95Lower,selected.CI95Upper-selected.MeanDifference,'-o');
    hold on;yline(0);xlabel('SNR (dB)');ylabel('Independent - shared excess MSE, CI95');grid on;
    export_figure(fig,output,'fig_r22d_test_difference_vs_snr');
    fig=new_figure();hold on;
    for index=1:3
        means=zeros(1,5);for velocity=1:5,means(velocity)=mean(channel_average(delta((velocity-1)*200+(1:200),:,primary_methods(index)),[]));end
        plot(1:5,means,'-o','Color',colors(index,:),'DisplayName',labels{index});
    end
    xticks(1:5);xticklabels(p.velocity_labels);ylabel('Regime-balanced excess MSE');grid on;legend('Location','best');
    export_figure(fig,output,'fig_r22d_test_excess_mse_vs_velocity');
    fig=new_figure();selected=paired(paired.Comparison=="AXIS_MINUS_"+methods(6).Name&startsWith(paired.Group,"SNR_"),:);
    errorbar(p.snr_db,selected.MeanDifference,selected.MeanDifference-selected.CI95Lower,selected.CI95Upper-selected.MeanDifference,'-o');
    hold on;yline(0);xlabel('SNR (dB)');ylabel('Independent - parameter-swapped excess MSE, CI95');grid on;
    export_figure(fig,output,'fig_r22d_swap_control');
    fig=new_figure();tiledlayout(1,2,'TileSpacing','compact');
    for filter={'h','g'}
        nexttile;stem(0:3,methods(4).Nu.(filter{1}),'DisplayName','Shared');hold on;
        stem(0:3,methods(5).Nu.(filter{1}),'DisplayName','Doppler');stem(0:3,methods(5).Tau.(filter{1}),'DisplayName','Delay');
        xlabel('FIR tap');ylabel('Coefficient');title(filter{1});grid on;legend('Location','best');
    end
    export_figure(fig,output,'fig_r22d_filter_nu_tau');
    fig=new_figure();tiledlayout(1,3,'TileSpacing','compact');
    for index=1:3,nexttile;imagesc(methods(primary_methods(index)).Mask);axis image;xlabel('Input stream');ylabel('Output stream');title(labels{index});end
    colormap(gray);export_figure(fig,output,'fig_r22d_support_maps');
    fig=new_figure();hold on;for index=1:3,boxchart(repmat(index,1000,1),frob_values(:,primary_methods(index)),'BoxFaceColor',colors(index,:));end
    xticks(1:3);xticklabels(labels);ylabel('Frozen-support Frobenius error');grid on;export_figure(fig,output,'fig_r22d_frobenius');
    fig=new_figure();hold on;for index=1:3,boxchart(repmat(index,1000,1),k99(:,primary_methods(index)),'BoxFaceColor',colors(index,:));end
    xticks(1:3);xticklabels(labels);ylabel('Descriptive TEST K99');grid on;export_figure(fig,output,'fig_r22d_k99');
end

function overlay(points,labels)
    for index=1:3,scatter(points(index,1),points(index,2),65,'filled','DisplayName',labels{index});end
end

function fig=new_figure()
    fig=figure('Visible','off','Color','w','Position',[100,100,1100,650]);
end

function export_figure(fig,output,name)
    exportgraphics(fig,fullfile(output,[name,'.pdf']),'ContentType','vector');
    exportgraphics(fig,fullfile(output,[name,'.png']),'Resolution',180);close(fig);
end

function flags=synthetic_check(p,output)
    smoke=fullfile(output,'synthetic_preflight');if ~isfolder(smoke),mkdir(smoke);end
    names=["D4_2D_PARTIAL_MALLAT","R2.2A_PUFB2D_STAR","R2.2C_CONFIRMED_STAR","R22D_SHARED_CONTROL","R22D_AXIS_SPECIFIC_STAR","R22D_AXIS_SWAPPED_CONTROL"];
    methods=repmat(r22d_model('transform',.78,.78),1,6);
    for method=1:6,methods(method).Name=names(method);methods(method).Mask=false(288);methods(method).Mask(1:p.K)=true;end
    methods(5)=r22d_model('transform',.7,.84);methods(5).Name=names(5);methods(5).Mask=methods(4).Mask;
    methods(6)=r22d_model('transform',.84,.7);methods(6).Name=names(6);methods(6).Mask=methods(4).Mask;
    rows=cell(48000,8);full_rows=cell(48000,5);frob_rows=cell(6000,5);sparse_rows=cell(6000,11);stability_rows=cell(30,6);
    cursor=0;channel_cursor=0;factors=[1,.8,.7,.6,.5,.55];
    for velocity=1:5
        for channel=1:200
            for method=1:6
                channel_cursor=channel_cursor+1;frob_rows(channel_cursor,:)={names(method),velocity,channel,p.K,.2*factors(method)};
                sparse_rows(channel_cursor,:)={names(method),velocity,channel,p.K,20000,25000,30000,.001,.8,.7,.85};
                for point=1:8
                    cursor=cursor+1;delta=.1*factors(method)+channel*method*1e-8;
                    rows(cursor,:)={names(method),velocity,channel,3*point,p.K,.1+delta,.1,delta};
                    full_rows(cursor,:)={names(method),velocity,channel,3*point,.1};
                end
            end
        end
        for method=1:6,stability_rows{(velocity-1)*6+method,1}=names(method);stability_rows((velocity-1)*6+method,2:6)={velocity,.7,.65,40,p.K};end
    end
    mse=cell2table(rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','RetainedEntries','MSESparse','MSEFull','DeltaMSE'});
    full=cell2table(full_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','FullMSE'});
    frob=cell2table(frob_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','RetainedEntries','FrobeniusError'});
    sparsity=cell2table(sparse_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','RetainedEntries','K95','K99','K999','J4','NormalizedEntropy','Hoyer','SupportCoverage'});
    stability=cell2table(stability_rows,'VariableNames',{'Method','VelocityIndex','WithinVelocityJaccard','CrossVelocityJaccard','DiagnosticChannelsPerVelocity','RetainedEntries'});
    alpha=(0:60)'*pi/61;coarse=table((1:61)',alpha,alpha,.5+(alpha-.78).^2,'VariableNames',{'CandidateID','AlphaNu','AlphaTau','Jrx'});
    writetable(coarse,fullfile(smoke,'r22d_shared_coarse_search.csv'));writetable(coarse,fullfile(smoke,'r22d_shared_refinement.csv'));
    [nu,tau]=ndgrid((0:16)*pi/17,(0:16)*pi/17);surface=table((1:289)',nu(:),tau(:),.4+(nu(:)-.7).^2+(tau(:)-.84).^2, ...
        'VariableNames',{'CandidateID','AlphaNu','AlphaTau','Jrx'});
    writetable(surface,fullfile(smoke,'r22d_axis_global_surface.csv'));writetable(surface,fullfile(smoke,'r22d_axis_refinement.csv'));
    selection=table((1:8)',[repmat("SHARED",3,1);repmat("AXIS",5,1)],repmat(.78,8,1),repmat(.78,8,1),(.5:.01:.57)', ...
        'VariableNames',{'CandidateID','Arm','AlphaNu','AlphaTau','ValJrx'});writetable(selection,fullfile(smoke,'r22d_val_selection.csv'));
    flags=r22d_finish(methods,p,mse,frob,sparsity,full,stability,smoke);
    assert(strcmp(flags.AXIS_SPECIFIC_MSE_BENEFIT,'SUPPORTED')&&flags.AXIS_SPECIALIZATION_SUPPORTED);
    r22a_io('json',fullfile(output,'r22d_preflight_report.json'),struct('PASS',true,'PhysicalChannelsGenerated',false,'Figures',12,'Comparisons',90));
end