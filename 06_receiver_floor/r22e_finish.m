function flags=r22e_finish(methods,p,output)
    if ischar(methods)&&strcmp(methods,'PLOT_TEST_MSE_ONLY')
        params=jsondecode(fileread(fullfile(output,'r22e_preregistration.json')));
        table=readtable(fullfile(output,'r22e_test_exact_mse.csv'),'TextType','string');
        exact=struct('Matrices',struct());
        for name=unique(table.Variant,'stable').'
            selected=sortrows(table(table.Variant==name,:),{'VelocityIndex','ChannelID','SNRdB'});
            exact.Matrices.(char(name))=reshape(selected.MSE,17,1000).';
        end
        destination=p;line_width=[];
        if isstruct(p),destination=p.OutputDirectory;line_width=p.LineWidth;end
        test_exact_mse_figure(params,exact,destination,line_width);return;
    end
    if ischar(methods)&&strcmp(methods,'PLOT_VELOCITY_ONLY')
        summary=readtable(fullfile(output,'r22e_ber_summary.csv'),'TextType','string');
        params=jsondecode(fileread(fullfile(output,'r22e_preregistration.json')));
        velocity_ber_figure(summary,params,p);return;
    end
    if ischar(methods)&&strcmp(methods,'PLOT_AGGREGATE_ONLY')
        summary=readtable(fullfile(output,'r22e_ber_summary.csv'),'TextType','string');
        aggregate_ber_figure(summary,p);return;
    end
    if ischar(methods)&&strcmp(methods,'SELFTEST'),flags=selftest(p,output);return;end
    exact=r22e_stats('MSE',p,output);ber=r22e_stats('BER',p,output);
    loaded=load(fullfile(output,'r22e_eta_freeze.mat'),'etas','common_eta');etas=loaded.etas;common_eta=loaded.common_eta;
    validation=r22e_math('read',fullfile(output,'r22e_val_confirmation.csv'));
    val=validation(validation.Method==methods(3).Name&validation.Region=="AGGREGATE",:);assert(height(val)==1);
    primary=exact.Comparisons(exact.Comparisons.Comparison=="AXIS_REG_MINUS_AXIS_ETA0"&exact.Comparisons.Region=="AGGREGATE",:);
    representation=exact.Comparisons(exact.Comparisons.Comparison=="AXIS_REG_MINUS_SHARED_REG"&exact.Comparisons.Region=="AGGREGATE",:);
    common=exact.Comparisons(exact.Comparisons.Comparison=="AXIS_COMMON_MINUS_SHARED_COMMON"&exact.Comparisons.Region=="AGGREGATE",:);
    numerical=jsondecode(fileread(fullfile(output,'r22e_numerical_controls.json')));
    WSE=mechanism_row(exact.Mechanism,'WSE_F');W=mechanism_row(exact.Mechanism,'W_F');omit=mechanism_row(exact.Mechanism,'T_OMIT');noise=mechanism_row(exact.Mechanism,'T_NOISE');
    mse_benefit='NOT_SUPPORTED';if primary.CI95Upper<0,mse_benefit='SUPPORTED';end
    mechanism=primary.CI95Upper<0&&WSE.CI95Upper<0&&omit.CI95Upper<0&&numerical.NUMERICAL_SOLVER_PASS;
    axis_shape=ber.Shapes(ber.Shapes.Group=="AGGREGATE"&ber.Shapes.Variant=="AXIS_REG",:);
    original_shape=ber.Shapes(ber.Shapes.Group=="AGGREGATE"&ber.Shapes.Variant=="AXIS_ETA0",:);
    high=ber.HighBER(ber.HighBER.Group=="AGGREGATE",:);
    corrected=axis_shape.HighSNRLogSlope<original_shape.HighSNRLogSlope&&high.CI95Upper<0&&primary.CI95Upper<0 ...
        &&WSE.CI95Upper<0&&numerical.NUMERICAL_SOLVER_PASS;
    robust=representation.CI95Upper<0&&common.CI95Upper<0;
    gain1=gain_row(ber.Gains,'AGGREGATE','SHARED_REG',1e-3);gain2=gain_row(ber.Gains,'AGGREGATE','SHARED_REG',1e-2);
    gain4=gain_row(ber.Gains,'AGGREGATE','SHARED_REG',1e-4);receiver_gain=gain_row(ber.Gains,'AGGREGATE','AXIS_ETA0',1e-3);
    high_gain=gain_row(ber.Gains,'HIGH_MOBILITY','SHARED_REG',1e-3);
    if etas(3)==0
        interpretation='NO_SUPPORTED_NONZERO_REGULARIZATION_FLOOR';
    elseif primary.CI95Upper>=0
        interpretation='REGULARIZATION_EFFECT_NOT_CONFIRMED_OUT_OF_SAMPLE';
    elseif mechanism&&corrected&&axis_shape.Shape=="CLEAR_DESCENT"
        interpretation='MISMATCH_AWARE_DIAGONAL_LOADING_CORRECTS_HIGH_SNR_SPARSE_RECEIVER_DEGRADATION';
    elseif (gain1.Support=="SUPPORTED"||receiver_gain.Support=="SUPPORTED")&&~mechanism
        interpretation='BER_IMPROVEMENT_WITHOUT_CONFIRMED_MISMATCH_MECHANISM';
    else
        interpretation='SECOND_ORDER_RECEIVER_GAIN_WITHOUT_RESOLVED_BER_GAIN';
    end
    flags=struct('M_FB',2,'NEW_TRANSFORM_TRAINED',false,'ALPHA_CHANGED',false,'FILTER_LENGTH_ESCALATED',false, ...
        'K_CHANGED',false,'SUPPORT_BUDGET_CHANGED',false,'RECEIVER_FAMILY_CHANGED',true, ...
        'RECEIVER_CHANGE','SINGLE_SCALE_NORMALIZED_DIAGONAL_LOADING_FLOOR','BER_USED_FOR_ETA_SELECTION',false, ...
        'TEST_USED_FOR_ETA_SELECTION',false,'HISTORICAL_BER_EXT_USED_FOR_SELECTION',false, ...
        'ETA_AXIS',etas(3),'ETA_SHARED',etas(2),'ETA_D4',etas(1),'ETA_R22A',etas(4),'ETA_COMMON',common_eta, ...
        'REGULARIZATION_VAL_SUPPORTED',val.CI95Upper<0,'TEST_MSE_REGULARIZATION_BENEFIT',mse_benefit, ...
        'MECHANISM_CORRECTION_CONFIRMED',mechanism,'HIGH_SNR_UPTURN_CORRECTED',corrected, ...
        'AXIS_REG_HIGH_SNR_SHAPE',char(axis_shape.Shape),'AXIS_ETA0_HIGH_SNR_SHAPE',char(original_shape.Shape), ...
        'BER_SNR_GAIN_1E2',char(gain2.Support),'BER_SNR_GAIN_1E3',char(gain1.Support),'BER_SNR_GAIN_1E4',char(gain4.Support), ...
        'AXIS_REGULARIZATION_GAIN_1E3',receiver_gain.MeasuredGainDb,'AXIS_VS_SHARED_GAIN_1E3',gain1.MeasuredGainDb, ...
        'HIGH_MOBILITY_BER_GAIN',char(high_gain.Support),'REPRESENTATION_ADVANTAGE_ROBUST_TO_REGULARIZATION',robust, ...
        'FULL_MSE_MONOTONIC_PASS',numerical.FULL_MSE_MONOTONIC_PASS,'NUMERICAL_SOLVER_PASS',numerical.NUMERICAL_SOLVER_PASS, ...
        'MANUSCRIPT_MODIFIED',false,'PRIMARY_INTERPRETATION',interpretation);
    r22a_io('json',fullfile(output,'r22e_decision_flags.json'),flags);
    train=r22e_math('read',fullfile(output,'r22e_train_eta_objective.csv'));baseline=train(train.Method==methods(3).Name&train.Eta==0,:);
    selected=train(train.Method==methods(3).Name&train.Eta==etas(3),:);
    receiver_gain2=gain_row(ber.Gains,'AGGREGATE','AXIS_ETA0',1e-2);receiver_gain4=gain_row(ber.Gains,'AGGREGATE','AXIS_ETA0',1e-4);
    lines=["# R2.2E Mismatch-Aware Regularized Sparse LMMSE";"1. M_FB:2.";"2. K:33178."; ...
        sprintf('3. alpha_nu:%.15g.',methods(3).AlphaNu);sprintf('4. alpha_tau:%.15g.',methods(3).AlphaTau); ...
        "5. TRAIN_RX:300 (60pervelocity).";"6. VAL_RX:200 (40pervelocity).";"7. TEST_RX:1000 (200pervelocity)."; ...
        "8. Frozen eta grid:"+string(mat2str(p.eta_grid)); ...
        sprintf('9. eta_axis_star:%.12g.',etas(3));sprintf('10. eta_shared_star:%.12g.',etas(2)); ...
        sprintf('11. eta_D4_star:%.12g.',etas(1));sprintf('12. eta_R22A_star:%.12g.',etas(4));sprintf('13. eta_common_star:%.12g.',common_eta); ...
        sprintf('14. TRAIN axis eta0 J_RX:%.12g.',baseline.TrainJrx);sprintf('15. TRAIN axis etaStar J_RX:%.12g.',selected.TrainJrx); ...
        sprintf('16. VAL axis etaStar-minus-eta0 MSE:%.12g.',val.MeanDifference);sprintf('17. VAL CI95:[%.12g,%.12g].',val.CI95Lower,val.CI95Upper); ...
        "18. REGULARIZATION_VAL_SUPPORTED="+string(val.CI95Upper<0);sprintf('19. TEST axis eta0 MSE:%.12g.',primary.ReferenceMean); ...
        sprintf('20. TEST axis etaStar MSE:%.12g.',primary.CandidateMean);sprintf('21. TEST relativeMSEreduction:%.9g%%.',100*(primary.ReferenceMean-primary.CandidateMean)/primary.ReferenceMean); ...
        sprintf('22. TEST HIGH30:36 WSE_F change:%.12g CI95[%.12g,%.12g].',WSE.MeanDifference,WSE.CI95Lower,WSE.CI95Upper); ...
        sprintf('23. TEST HIGH W_F change:%.12g CI95[%.12g,%.12g].',W.MeanDifference,W.CI95Lower,W.CI95Upper); ...
        sprintf('24. TEST HIGH T_OMIT change:%.12g CI95[%.12g,%.12g].',omit.MeanDifference,omit.CI95Lower,omit.CI95Upper); ...
        sprintf('25. TEST HIGH T_NOISE change:%.12g CI95[%.12g,%.12g].',noise.MeanDifference,noise.CI95Lower,noise.CI95Upper); ...
        "26. MECHANISM_CORRECTION_CONFIRMED="+string(mechanism); ...
        sprintf('27. AXIS_ETA0 minimumBER:%.12g.',original_shape.MinimumBER);sprintf('28. AXIS_ETA0 minimumBER SNR:%d.',original_shape.MinimumBERSNR); ...
        sprintf('29. AXIS_REG minimumBER:%.12g.',axis_shape.MinimumBER);sprintf('30. AXIS_REG minimumBER SNR:%d.',axis_shape.MinimumBERSNR); ...
        sprintf('31. AXIS_ETA0 highSNRslope:%.12g decades/dB.',original_shape.HighSNRLogSlope); ...
        sprintf('32. AXIS_REG highSNRslope:%.12g decades/dB.',axis_shape.HighSNRLogSlope); ...
        "33. AXIS_ETA0 shape="+original_shape.Shape;"34. AXIS_REG shape="+axis_shape.Shape;"35. HIGH_SNR_UPTURN_CORRECTED="+string(corrected); ...
        sprintf('36. Axisregularization gain1e-2:%.12g dB (%s).',receiver_gain2.MeasuredGainDb,receiver_gain2.Support); ...
        sprintf('37. Axisregularization gain1e-3:%.12g dB (%s).',receiver_gain.MeasuredGainDb,receiver_gain.Support); ...
        sprintf('38. CI95:[%.12g,%.12g]; resolvedbootstrap%d/%d.',receiver_gain.CI95LowerDb,receiver_gain.CI95UpperDb,receiver_gain.ResolvedResamples,receiver_gain.Resamples); ...
        sprintf('39. Axisregularization gain1e-4:%.12g dB (%s).',receiver_gain4.MeasuredGainDb,receiver_gain4.Support); ...
        sprintf('40. AxisvsShared REG gain1e-3:%.12g dB (%s).',gain1.MeasuredGainDb,gain1.Support); ...
        sprintf('41. CI95:[%.12g,%.12g]; resolvedbootstrap%d/%d.',gain1.CI95LowerDb,gain1.CI95UpperDb,gain1.ResolvedResamples,gain1.Resamples)];
    for velocity=1:5
        row=gain_row(ber.Gains,sprintf('VELOCITY_%d',velocity),'SHARED_REG',1e-3);
        lines(end+1)=sprintf('%d. %s axis-vs-sharedREG gain1e-3:%.12g dB CI95[%.12g,%.12g] (%s).',41+velocity,p.velocity_labels{velocity}, ...
            row.MeasuredGainDb,row.CI95LowerDb,row.CI95UpperDb,row.Support); %#ok<AGROW>
    end
    lines=[lines;"47. HIGH_MOBILITY_BER_GAIN="+high_gain.Support; ...
        sprintf('48. Commoneta axis-minus-sharedMSE:%.12g CI95[%.12g,%.12g]; method-specific:%.12g CI95[%.12g,%.12g].',common.MeanDifference,common.CI95Lower,common.CI95Upper,representation.MeanDifference,representation.CI95Lower,representation.CI95Upper); ...
        "49. REPRESENTATION_ADVANTAGE_ROBUST_TO_REGULARIZATION="+string(robust); ...
        "50. FULL_MSE_MONOTONIC_PASS="+string(numerical.FULL_MSE_MONOTONIC_PASS);"51. NUMERICAL_SOLVER_PASS="+string(numerical.NUMERICAL_SOLVER_PASS); ...
        "52. BER_USED_FOR_ETA_SELECTION=FALSE";"53. TEST_USED_FOR_ETA_SELECTION=FALSE";"54. NEW_TRANSFORM_TRAINED=FALSE"; ...
        "55. SUPPORT_BUDGET_CHANGED=FALSE";"56. M_FB=2";"57. MANUSCRIPT_MODIFIED=FALSE";"58. PRIMARY_INTERPRETATION="+string(interpretation); ...
        "";"## Protocol And Claim Limits"; ...
        "Section7explicitlyrequiresNEW TRAIN_RX commonmasks using historicalranking withK33178. These masks are newlyconstructedfrom300unseenTRAINchannels and frozenBEFORE eta scoring; not the previouslyhistoricalmask. NoVAL/TEST masklearning; representations/alphas unchanged. Supportbudgetunchanged doesnot meanmaskunchangedversushistory."; ...
        "OnlyONEscalar eta perrepresentation; lambda=sigma2+eta*norm(Hs,F)^2/n operationallyusesHs,sigma2,etaonly. H/E only offlineexactobjective/diagnostics. The grid andsplits frozenbeforegeneration,eta selectedTRAINonly,VALconfirmationwithoutretuning."; ...
        "NohistoricalBEREXTchannels/scoreswereusedinselection. TESTopenedonlyafterallfourmethodetas/commoneta,supports,sources,statistics,BERstoppingandfigures frozen. BERusespairednew communicationrealizationson1000TEST_RXchannels."; ...
        "BERcounts pooledwithinvelocity,aggregateequallyweightedvelocities. Physicalchannelbootstrap10000pairedreplicates; targetcrossingsadjacentpositivemeasuredpointswith>=500aggregateerrors,linearSNRvslogBER,noextrapolation/splines.>=95%pairedresamplesmustresolve;conditionalCIandresolvedcountdisclosed. MissingcrossingsNaN/null arenotzero gains."; ...
        "ShapeOLS30:36positiveBERrequires>=3points; CLEAR_DESCENT<=-.03,APPROX_FLAT(-.03,.03),UPTURN>=.03. APPROX_FLATisnotasymptoticerrorfloor. Correctionflagalso requirespairedhighBER,MSE,WSEreductionandnumericalcontrols; strong wording requiresCLEAR_DESCENT."; ...
        "Controls:eta0reproduceshistoricalreceiverfamilyonnewmask/bank; fullchannelMSEwithouteta floor; commonetacontrolisolatesrepresentationalgainfrommethod-specifictuning. No claimofoptimalreceiver/optimalregularization/universalBERimprovement orcompleteeliminationofmismatch."; ...
        sprintf('Numericaldirect-risk maximumdifference%.12g; backwardmaximum%.12g.',numerical.MaxDirectSolveMSEError,numerical.MaxBackwardResidual)];
    r22a_io('text',fullfile(output,'r22e_final_report.md'),join(lines,newline));r22a_io('text',fullfile(output,'r22e_executive_summary.md'),join(lines,newline));
    figures(methods,p,exact,ber,train,validation,output);
end

function aggregate_ber_figure(summary,output)
    % Legend order; drawing order puts the primary receiver on top.
    variants=["AXIS_REG","SHARED_REG","D4_REG","R22A_REG","AXIS_ETA0"];
    draw_order=[3,5,2,4,1];
    % Marker phase (SNR mod 4) differs per curve so coincident curves stay distinguishable.
    marker_indices=struct('AXIS_REG',1:4:37,'AXIS_ETA0',2:4:37,'SHARED_REG',3:4:37,'R22A_REG',4:4:37,'D4_REG',1:4:37);
    [fig,ax]=paper_figure(7.16,3.80);
    yline(ax,1e-3,'--','Color',[0.55 0.55 0.55],'LineWidth',0.6,'HandleVisibility','off');
    curve_handles=gobjects(1,numel(variants));labels=cell(1,numel(variants));curves=struct([]);lower=inf;upper=0;
    for index=draw_order
        T=summary(summary.Group=="AGGREGATE"&summary.Variant==variants(index)&summary.SNRdB>=0,:);
        T=sortrows(T,'SNRdB');values=T.BER;values(values<=0)=NaN;
        style=receiver_style(variants(index));labels{index}=style.Label;
        bars=errorbar(ax,T.SNRdB,values,max(0,T.BER-T.Lower95),max(0,T.Upper95-T.BER), ...
            'LineStyle','none','Marker','none','Color',style.Color,'LineWidth',0.6, ...
            'CapSize',2.5,'HandleVisibility','off');
        curve_handles(index)=plot(ax,T.SNRdB,values,style.Arguments{:}, ...
            'MarkerIndices',marker_indices.(char(variants(index))));
        curves=[curves,plotted_curve(variants(index),style.Label,curve_handles(index),bars)]; %#ok<AGROW>
        lower=min([lower,bars.YData-bars.YNegativeDelta]);upper=max([upper,bars.YData+bars.YPositiveDelta]);
    end
    % Limits from the stored values and CIs, with ~0.1-decade margins.
    margin=0.04*log10(upper/lower);
    set(ax,'YScale','log','XLim',[0,36],'XTick',0:4:36,'YLim',[lower,upper].*10.^[-margin,margin], ...
        'YTick',[1e-3,1e-2,1e-1],'YMinorTick','on');
    paper_axes(ax);xlabel(ax,'SNR (dB)');ylabel(ax,'Measured BER');
    paper_legend(ax,curve_handles,labels);
    paper_export(fig,output,'fig_r22e_ber_aggregate',curves);
end

function velocity_ber_figure(summary,p,output)
    probe=new_figure();layout=tiledlayout(probe,2,3,'TileSpacing','compact');
    probe_axes=gobjects(1,5);
    for velocity=1:5
        probe_axes(velocity)=nexttile(layout);
        velocity_ber_panel(probe_axes(velocity),summary,p,velocity);
    end
    drawnow;
    positions=vertcat(probe_axes.Position);
    close(probe);

    column_spacing=positions(2,1)-positions(1,1);
    positions(4,1)=positions(1,1)+column_spacing/2;
    positions(5,1)=positions(2,1)+column_spacing/2;

    fig=new_figure();panel_axes=gobjects(1,5);
    for velocity=1:5
        panel_axes(velocity)=axes('Parent',fig,'Position',positions(velocity,:));
        velocity_ber_panel(panel_axes(velocity),summary,p,velocity);
    end
    drawnow;
    final_positions=vertcat(panel_axes.Position);
    assert(max(abs(final_positions(1:3,:)-positions(1:3,:)),[],'all')<1e-10);
    assert(max(abs(final_positions(4:5,3:4)-positions(1,3:4)),[],'all')<1e-10);
    assert(abs(mean(final_positions(4:5,1)+final_positions(4:5,3)/2)-0.5)<1e-10);
    export_figure(fig,output,'fig_r22e_ber_by_velocity');
end

function velocity_ber_panel(ax,summary,p,velocity)
    axes(ax);hold on;
    for name=["AXIS_REG","SHARED_REG","AXIS_ETA0"]
        plot_curve(summary,sprintf('VELOCITY_%d',velocity),name,0);
    end
    axes_style();title(p.velocity_labels{velocity});
    if velocity==1,legend('Interpreter','none','Location','best');end
end

function row=mechanism_row(T,field)
    row=T(T.Metric==string(field),:);assert(height(row)==1);
end

function row=gain_row(T,group,reference,target)
    row=T(T.Group==string(group)&T.Comparison=="AXIS_REG_VS_"+string(reference)&T.TargetBER==target,:);assert(height(row)==1);
end

function figures(methods,p,exact,ber,train,validation,output)
    fig=new_figure();hold on;
    for method=1:4,T=train(train.Method==methods(method).Name,:);plot(1:14,T.TrainJrx,'-o','DisplayName',figure_label(methods(method).Name));end
    xticks(1:14);xticklabels(compose('%.0e',p.eta_grid));xtickangle(45);ylabel('Equal-SNR TRAIN true MSE');xlabel('Frozen eta');grid on;legend('Interpreter','none','Location','best');
    export_figure(fig,output,'fig_r22e_train_eta_landscape');
    fig=new_figure();T=validation(validation.Region=="AGGREGATE",:);bar([T.Eta0MSE,T.EtaStarMSE]);xticks(1:4);xticklabels(arrayfun(@figure_label,T.Method,'UniformOutput',false));ylabel('VAL exact MSE');legend('eta=0','TRAIN selected eta');grid on;
    export_figure(fig,output,'fig_r22e_val_eta_confirmation');
    test_exact_mse_figure(p,exact,output,1.4);
    diagnostic=r22e_math('read',fullfile(output,'r22e_test_mismatch_decomposition.csv'));
    fig=new_figure();tiledlayout(1,3,'TileSpacing','compact');
    for field={'T_OMIT','T_NOISE','MismatchPenalty'}
        nexttile;hold on;for name=["AXIS_ETA0","AXIS_REG"],T=diagnostic(diagnostic.Variant==name,:);group=groupsummary(T,'SNRdB','mean',field{1});plot(group.SNRdB,group.(['mean_',field{1}]),'-o','DisplayName',figure_label(name));end
        xlabel('SNR (dB)');title(field{1},'Interpreter','none');grid on;legend('Interpreter','none');
    end
    export_figure(fig,output,'fig_r22e_mismatch_terms');
    fig=new_figure();tiledlayout(1,2,'TileSpacing','compact');
    for field={'W_F','WSE_F'}
        nexttile;hold on;for name=["AXIS_ETA0","AXIS_REG"],T=diagnostic(diagnostic.Variant==name,:);group=groupsummary(T,'SNRdB','mean',field{1});plot(group.SNRdB,group.(['mean_',field{1}]),'-o','DisplayName',figure_label(name));end
        xlabel('SNR (dB)');title(field{1},'Interpreter','none');grid on;legend('Interpreter','none');
    end
    export_figure(fig,output,'fig_r22e_receiver_gain');
    aggregate_ber_figure(ber.Summary,output);
    velocity_ber_figure(ber.Summary,p,output);
    fig=new_figure();hold on;for name=["AXIS_REG","AXIS_ETA0"],plot_curve(ber.Summary,'AGGREGATE',name,0);end
    axes_style();legend('Interpreter','none','Location','best');export_figure(fig,output,'fig_r22e_axis_regularization_control');
    fig=new_figure();tiledlayout(1,2,'TileSpacing','compact');groups=["AGGREGATE","VELOCITY_1","VELOCITY_2","VELOCITY_3","VELOCITY_4","VELOCITY_5","HIGH_MOBILITY"];
    for reference=["SHARED_REG","AXIS_ETA0"]
        nexttile;hold on;
        for index=1:7
            T=gain_row(ber.Gains,groups(index),reference,1e-3);
            if T.Support~="NOT_RESOLVED",errorbar(index,T.MeasuredGainDb,T.MeasuredGainDb-T.CI95LowerDb,T.CI95UpperDb-T.MeasuredGainDb,'o');
            else,text(index,0,'NOT RESOLVED','Rotation',45);end
        end
        yline(0);xlim([.5,7.5]);xticks(1:7);xticklabels({'All','50','100','500','Mach1','Mach5','High'});ylabel('SNR gain at BER1e-3 (dB)');title([figure_label('AXIS_REG'),' vs. ',figure_label(reference)],'Interpreter','none');grid on;
    end
    export_figure(fig,output,'fig_r22e_snr_gain_1e3');
    fig=new_figure();hold on;for name=["AXIS_REG","AXIS_ETA0"],plot_curve(ber.Summary,'AGGREGATE',name,30);end
    axes_style();legend('Interpreter','none','Location','best');export_figure(fig,output,'fig_r22e_highsnr_shape');
end

function test_exact_mse_figure(p,exact,output,~)
    % Legend order; drawing order puts the primary receiver on top.
    variants=["AXIS_REG","SHARED_REG","AXIS_ETA0","D4_REG","R22A_REG"];
    draw_order=[4,3,2,5,1];
    marker_indices=struct('AXIS_REG',1:2:17,'SHARED_REG',1:2:17,'AXIS_ETA0',1:2:17,'D4_REG',2:2:17,'R22A_REG',2:2:17);
    [fig,ax]=paper_figure(7.16,3.70);
    curve_handles=gobjects(1,numel(variants));labels=cell(1,numel(variants));curves=struct([]);
    for index=draw_order
        style=receiver_style(variants(index));labels{index}=style.Label;
        curve_handles(index)=plot(ax,p.development_snr,mean(exact.Matrices.(char(variants(index))),1), ...
            style.Arguments{:},'MarkerIndices',marker_indices.(char(variants(index))));
        curves=[curves,plotted_curve(variants(index),style.Label,curve_handles(index),[])]; %#ok<AGROW>
    end
    % Limits from the stored curves with a 6% margin, rounded outward to 0.005.
    values=[curves.Y];span=max(values)-min(values);
    limits=[floor((min(values)-0.06*span)/0.005),ceil((max(values)+0.06*span)/0.005)]*0.005;
    set(ax,'XLim',[20,36],'XTick',20:2:36,'YLim',limits);ax.YAxis.TickLabelFormat='%.2f';
    paper_axes(ax);xlabel(ax,'SNR (dB)');ylabel(ax,'True-channel MSE');
    paper_legend(ax,curve_handles,labels);
    paper_export(fig,output,'fig_r22e_test_exact_mse',curves);
end

function style=receiver_style(name)
    % One method -> color/line/marker identity shared by Figs. 4 and 6; redundant
    % line style and marker keep the encoding readable in grayscale.
    switch char(name)
        case 'AXIS_REG',  label='Independent, \eta = 0.01';color=[15,77,146]/255; line='-'; marker='o';face='filled';width=2.0;size=6;
        case 'SHARED_REG',label='Shared, \eta = 0.01';     color=[182,67,66]/255; line='-'; marker='s';face='filled';width=1.5;size=6.5;
        case 'AXIS_ETA0', label='Independent, \eta = 0';   color=[154,77,142]/255;line='--';marker='^';face='open';  width=1.5;size=6.5;
        case 'D4_REG',    label='D4, \eta = 0.03';         color=[77,77,77]/255;  line='-.';marker='d';face='open';  width=1.3;size=6.5;
        case 'R22A_REG',  label='Compressibility-selected shared';color=[66,148,158]/255;line=':';marker='v';face='open';width=1.8;size=6.5;
        otherwise,error('No receiver style for %s.',name);
    end
    fill=[1,1,1];if strcmp(face,'filled'),fill=color;end
    style=struct('Label',label,'Color',color,'Arguments',{{'LineStyle',line,'Marker',marker,'Color',color, ...
        'LineWidth',width,'MarkerSize',size,'MarkerFaceColor',fill,'MarkerEdgeColor',color}});
end

function [fig,ax]=paper_figure(width,height)
    % Source geometry is the final PDF page size; text is sized for ~0.72\textwidth.
    fig=figure('Visible','off','Color','w','Units','inches','Position',[0.5,0.5,width,height], ...
        'PaperUnits','inches','PaperSize',[width,height],'PaperPosition',[0,0,width,height], ...
        'PaperPositionMode','manual','InvertHardcopy','off');
    ax=axes(fig);hold(ax,'on');
end

function paper_axes(ax)
    set(ax,'FontName','Helvetica','FontSize',11,'LineWidth',0.75,'Box','off','TickDir','out', ...
        'TickLength',[0.007,0.007],'XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off', ...
        'GridLineStyle','-','GridColor',[0.15,0.15,0.15],'GridAlpha',0.12,'Layer','bottom', ...
        'XColor',[0.15,0.15,0.15],'YColor',[0.15,0.15,0.15],'TickLabelInterpreter','tex');
    ax.XLabel.FontSize=12;ax.YLabel.FontSize=12;
end

function paper_legend(ax,handles,labels)
    % Column-major fill keeps the loaded Independent/Shared pair in the first column.
    legend_handle=legend(ax,handles,labels,'Interpreter','tex','Orientation','vertical','NumColumns',3, ...
        'FontName','Helvetica','FontSize',11,'Box','off','Location','northoutside');
    legend_handle.ItemTokenSize=[34,18];
    % Legend band across the top, axes filling the rest with a small gap (inches).
    drawnow;set([ax,legend_handle],'Units','inches');page=ax.Parent.Position(3:4);
    box=legend_handle.Position;box(1:2)=[(page(1)-box(3))/2,page(2)-box(4)-0.02];
    legend_handle.Position=box;inset=ax.TightInset;
    ax.Position=[inset(1)+0.04,inset(2)+0.07,page(1)-inset(1)-0.04-0.16,box(2)-0.10-inset(2)-0.07];
end

function curve=plotted_curve(name,label,line,bars)
    % Values read back from the rendered graphics objects, for data-identity checks.
    curve=struct('Variant',string(name),'Label',string(label),'X',line.XData,'Y',line.YData, ...
        'YNegativeDelta',[],'YPositiveDelta',[]);
    if ~isempty(bars)
        assert(isequaln(bars.XData,line.XData)&&isequaln(bars.YData,line.YData));
        curve.YNegativeDelta=bars.YNegativeDelta;curve.YPositiveDelta=bars.YPositiveDelta;
    end
end

function paper_export(fig,output,name,curves)
    drawnow;
    print(fig,fullfile(output,[name,'.pdf']),'-dpdf','-vector');
    print(fig,fullfile(output,[name,'.png']),'-dpng','-r600');
    save(fullfile(output,[name,'_plotted_data.mat']),'curves');
    close(fig);
end

function plot_curve(summary,group,name,minimum)
    T=summary(summary.Group==string(group)&summary.Variant==name&summary.SNRdB>=minimum,:);values=T.BER;values(values<=0)=NaN;
    errorbar(T.SNRdB,values,max(0,T.BER-T.Lower95),max(0,T.Upper95-T.BER),'-o','DisplayName',figure_label(name));
    zero=T.BER==0;if any(zero),scatter(T.SNRdB(zero),T.Upper95(zero),25,'v','DisplayName',[figure_label(name),' zero-error upperbound']);end
end

function label=figure_label(name)
    switch char(name)
        case {'AXIS_REG'},label='Independent, eta=0.01';
        case {'AXIS_ETA0'},label='Independent, eta=0';
        case {'SHARED_REG'},label='Shared, eta=0.01';
        case {'D4_REG'},label='D4, eta=0.03';
        case {'R22A_REG','R2.2A_PUFB2D_STAR'},label='Compressibility-trained shared';
        case {'R2.2D_AXIS_RX'},label='Independent';
        case {'R2.2C_SHARED_RX'},label='Shared';
        case {'D4_2D_PARTIAL_MALLAT'},label='D4';
        otherwise,label=char(name);
    end
end

function axes_style()
    set(gca,'YScale','log');grid on;xlabel('SNR (dB)');ylabel('Measured BER, confidence interval');
end

function fig=new_figure()
    fig=figure('Visible','off','Color','w','Position',[100,100,1100,650]);
end

function export_figure(fig,output,name)
    exportgraphics(fig,fullfile(output,[name,'.pdf']),'ContentType','vector');exportgraphics(fig,fullfile(output,[name,'.png']),'Resolution',180);close(fig);
end

function flags=selftest(p,output)
    smoke=fullfile(output,'synthetic_preflight');if ~isfolder(smoke),mkdir(smoke);end
    loaded=load(fullfile(output,'r22e_frozen_transforms.mat'),'methods');methods=loaded.methods;
    etas=[.03,.01,.01,.01];common_eta=.01;save(fullfile(smoke,'r22e_eta_freeze.mat'),'etas','common_eta');
    risk_rows=cell(56,3);cursor=0;
    for method=1:4,for eta=p.eta_grid,cursor=cursor+1;risk_rows(cursor,:)={methods(method).Name,eta,.1+(eta-.01)^2};end;end
    writetable(cell2table(risk_rows,'VariableNames',{'Method','Eta','TrainJrx'}),fullfile(smoke,'r22e_train_eta_objective.csv'));
    val_rows=cell(16,7);cursor=0;
    for method=1:4,for region=["AGGREGATE","LOW","MID","HIGH"],cursor=cursor+1;val_rows(cursor,:)={methods(method).Name,region,.1,.05,-.05,-.06,-.04};end;end
    writetable(cell2table(val_rows,'VariableNames',{'Method','Region','Eta0MSE','EtaStarMSE','MeanDifference','CI95Lower','CI95Upper'}),fullfile(smoke,'r22e_val_confirmation.csv'));
    variants=["AXIS_ETA0","AXIS_REG","SHARED_ETA0","SHARED_REG","D4_REG","R22A_REG","D4_COMMON","SHARED_COMMON","AXIS_COMMON","R22A_COMMON"];
    rows=cell(170000,7);diag_rows=cell(34000,15);cursor=0;diag_cursor=0;
    for channel=1:1000
        velocity=ceil(channel/200);local=mod(channel-1,200)+1;
        for variant=1:10
            for point=1:17
                cursor=cursor+1;value=.1;if ismember(variant,[2,9]),value=.05;end
                rows(cursor,:)={variants(variant),velocity,local,19+point,.01,33178,value+channel*1e-8};
                if variant<=2
                    diag_cursor=diag_cursor+1;scale=1;if variant==2,scale=.5;end
                    diag_rows(diag_cursor,:)={variants(variant),velocity,local,19+point,value,.01,.05*scale,0,.01*scale,10*scale,5*scale,100,1,.9,.05*scale};
                end
            end
        end
    end
    writetable(cell2table(rows,'VariableNames',{'Variant','VelocityIndex','ChannelID','SNRdB','Eta','RetainedEntries','MSE'}),fullfile(smoke,'r22e_test_exact_mse.csv'));
    writetable(cell2table(diag_rows,'VariableNames',{'Variant','VelocityIndex','ChannelID','SNRdB','MSE','T_REG','T_OMIT','T_CROSS','T_NOISE','W_F','WSE_F','CondA','LambdaEff','Qs','MismatchPenalty'}),fullfile(smoke,'r22e_test_mismatch_decomposition.csv'));
    r22a_io('json',fullfile(smoke,'r22e_numerical_controls.json'),struct('NUMERICAL_SOLVER_PASS',true,'FULL_MSE_MONOTONIC_PASS',true,'MaxDirectSolveMSEError',0,'MaxBackwardResidual',0));
    rows=cell(185000,9);cursor=0;names=["AXIS_REG","SHARED_REG","D4_REG","R22A_REG","AXIS_ETA0"];shifts=[2,0,-1,0,0];
    for velocity=1:5
        for channel=1:200
            for variant=1:5
                for snr=0:36
                    rate=.1*10^(-(snr+shifts(variant))/10);
                    if variant==5&&snr>24,rate=.1*10^(-24/10)*10^(.05*(snr-24));end
                    cursor=cursor+1;rows(cursor,:)={names(variant),velocity,channel,snr,round(1e6*rate),1e6,1736,33178,.01};
                end
            end
        end
    end
    writetable(cell2table(rows,'VariableNames',{'Variant','VelocityIndex','ChannelID','SNRdB','Errors','Bits','Frames','RetainedEntries','Eta'}),fullfile(smoke,'r22e_ber_per_channel.csv'));
    flags=r22e_finish(methods,p,smoke);assert(flags.MECHANISM_CORRECTION_CONFIRMED&&flags.HIGH_SNR_UPTURN_CORRECTED ...
        &&abs(flags.AXIS_VS_SHARED_GAIN_1E3-2)<.02&&flags.REPRESENTATION_ADVANTAGE_ROBUST_TO_REGULARIZATION);
    r22a_io('json',fullfile(output,'r22e_report_preflight.json'),struct('PASS',true,'PhysicalChannelsGenerated',false,'Known2dBGainReproduced',true,'Figures',10));
end