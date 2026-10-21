function result=r29_stage_stats(paths,p)
% TEST statistics, predeclared gates (r29_candidate_freeze.mat: gates) and M4-A/B/C/D classification.
% Equal-velocity means; paired channel bootstrap within velocity (10000 replicates, seed 292200000), percentile 95% CIs.
    out=fullfile(paths.d,'data');a=load(fullfile(out,'r29d_analytic.mat'));s=load(fullfile(out,'r29d_support.mat'));
    g=load(fullfile(paths.c,'config','r29_candidate_freeze.mat'),'gates');g=g.gates;
    names=string(a.names);rho=p.rho_test;nK=numel(rho);count=size(a.terms,1);vel=repelem((1:5)',count/5);
    agg=@(x)(mean(x(:,:,:,1:3),4)+mean(x(:,:,:,4:6),4)+mean(x(:,:,:,7:8),4))/3;     % arm R regime-balanced (as R2.6)
    full=reshape(a.full,count,1,1,[]);
    M.DeltaM=agg(a.terms(:,:,:,:,5,1)-full);M.DeltaMOracle=agg(a.terms(:,:,:,:,5,2)-full);M.RDeltaM=M.DeltaM-M.DeltaMOracle;
    M.TOmit=agg(a.terms(:,:,:,:,2,1));M.TCross=agg(a.terms(:,:,:,:,3,1));M.TNoise=agg(a.terms(:,:,:,:,4,1));M.TReg=agg(a.terms(:,:,:,:,1,1));
    M.OmittedOffline=a.scalar(:,:,:,1,1);M.OmittedOracle=a.scalar(:,:,:,1,2);M.CapturedOffline=a.scalar(:,:,:,4,1);M.CapturedOracle=a.scalar(:,:,:,4,2);
    M.RE=M.OmittedOffline-M.OmittedOracle;M.Overlap=s.overlap;M.EtaChannel=s.eta;M.ReceiverNorm12dB=a.scalar(:,:,:,2,1);
    amp=M.TOmit./M.OmittedOffline;amp(M.OmittedOffline<=1e-15)=NaN;M.MismatchAmplification=amp;
    rng(292200000,'twister');B=10000;W=zeros(B,count);for v=1:5,idx=find(vel==v);d=randi(numel(idx),B,numel(idx));
        for b=1:B,W(b,idx)=accumarray(d(b,:)',1,[numel(idx),1])';end,end
    W=W/(count/5)/5;w0=ones(1,count)/count;   % equal-velocity weights (balanced design: 200 per velocity)
    % ---------- curves
    rows={};fields=fieldnames(M);
    for r=1:numel(names),for k=1:nK
        row={names(r),string(a.groups{r}),rho(k)};for f=fields',x=M.(f{1})(:,r,k);row{end+1}=mean(x,'omitnan');end %#ok<AGROW>
        row=[row,{mean(M.CapturedOffline(:,r,k))/mean(M.CapturedOracle(:,r,k)),median(s.eta(:,r,k)),prctile(s.eta(:,r,k),10),prctile(s.eta(:,r,k),90), ...
            double(s.Hs(r,k)),mean(s.jac(:,r,k)),mean(s.jacv(:,:,r,k),'all')}]; %#ok<AGROW>
        rows(end+1,:)=row; %#ok<AGROW>
    end,end
    C=cell2table(rows,'VariableNames',[{'Representation','Group','Rho'},fields',{'EtaPooled','EtaMedian','EtaP10','EtaP90','HSupport','JaccardOverall','JaccardWithinVelocity'}]);
    writetable(C,fullfile(out,'r29d_curves.csv'));
    % ---------- comparisons vs best comparator (per rho, per metric)
    comp=find(ismember(names,string(g.comparator_set)));tested=find(ismember(names,string(g.tested_architectures)));
    metrics={'DeltaM','OmittedOffline','OmittedOracle','RE','EtaPooled','JaccardOverall'};rows={};
    for m=metrics
        for k=1:nK
            [ctrl,cv]=best_control(m{1},M,s,comp,k,w0);
            for t=[tested(:)',find(names=="HYBRID-M4-DYADIC")]
                [est,lo,hi]=relative(m{1},M,s,t,ctrl,k,W,w0);
                rows(end+1,:)={string(m{1}),rho(k),names(t),names(ctrl),cv,est,lo,hi,lo>0}; %#ok<AGROW>
            end
        end
    end
    R=cell2table(rows,'VariableNames',{'Metric','Rho','Representation','BestComparator','ComparatorMean','RelativeImprovement','Lo','Hi','Material'});
    R.Material=R.Material&R.RelativeImprovement>=0.10&~(R.Metric=="EtaPooled"|R.Metric=="JaccardOverall");
    R.Material(R.Metric=="EtaPooled"|R.Metric=="JaccardOverall")=R.Lo(R.Metric=="EtaPooled"|R.Metric=="JaccardOverall")>0;   % increase with CI > 0
    writetable(R,fullfile(out,'r29d_vs_best_comparator.csv'));
    % ---------- predeclared gates and classification (tested architectures only; HYBRID diagnostic)
    grid=g.practical_grid;cls=struct();berA=false;berB=false;
    for t=tested(:)'
        sel=@(metric)arrayfun(@(r)any(R.Material(R.Metric==metric&R.Representation==names(t)&abs(R.Rho-r)<1e-9)),grid);
        dm=sel("DeltaM");off=sel("OmittedOffline");or=sel("OmittedOracle");re=sel("RE");eta=sel("EtaPooled");jac=sel("JaccardOverall");
        re25=arrayfun(@(r)any(R.RelativeImprovement(R.Metric=="RE"&R.Representation==names(t)&abs(R.Rho-r)<1e-9)>=0.25&R.Lo(R.Metric=="RE"&R.Representation==names(t)&abs(R.Rho-r)<1e-9)>0),grid);
        stable=arrayfun(@(r)velocity_sign(M,comp,t,find(abs(rho-r)<1e-9),w0),grid);
        A=runs(dm)>=2;Bg=runs(re25&eta&stable)>=2;berA=berA||A;berB=berB||Bg;
        practical=off|eta|re;
        if ~any(dm|off|or|re|eta|jac),c='M4-A';elseif runs(practical)>=3,c='M4-C';elseif runs(or)>=2&&runs(practical)<2,c='M4-B';else,c='M4-D';end
        cls.(matlab.lang.makeValidName(names(t)))=struct('class',c,'DeltaM',dm,'Offline',off,'Oracle',or,'RE',re,'Eta',eta,'Jaccard',jac,'RE25',re25,'VelocityStable',stable,'BER_A',A,'BER_B',Bg);
    end
    % ---------- crossover structure (best mean DeltaM per rho; Identity vs PUFB, PUFB vs OTFS)
    cross=table(rho(:),'VariableNames',{'Rho'});
    for k=1:nK,[~,b]=min(mean(M.DeltaM(:,[comp(:);tested(:)],k),1));all_=[comp(:);tested(:)];cross.Best(k)=names(all_(b));
        [~,b]=min(mean(M.DeltaM(:,comp,k),1));cross.BestComparator(k)=names(comp(b));end
    writetable(cross,fullfile(out,'r29d_best_by_rho.csv'));
    % ---------- causal pairs (sec. 22): difference (second - first) in each metric with CI; negative = second lower
    pairs={'A_support8','DYADIC-W1-D3','M4-BASIC-D1';'A_support16','DYADIC-W1-D7','M4-BASIC-(3,3)';'B_structured_vs_basic','M4-BASIC-(3,3)','M4-INDEPENDENT-STRUCT-(1,1)'; ...
        'C_axis_capacity_depth','M4-INDEPENDENT-STRUCT-(1,1)','M4-INDEPENDENT-STRUCT-(1,2)';'D_terminal_resolution','M4-INDEPENDENT-STRUCT-(1,1)','HYBRID-M4-DYADIC'; ...
        'E_axis_sharing','M4-INDEPENDENT-STRUCT-(1,1)','M4-SHARED-STRUCT-(1,1)';'X_current_vs_reoptimized_D1','INDEPENDENT_PUFB','DYADIC-W1-D1'; ...
        'X_D1_vs_M4_struct','DYADIC-W1-D1','M4-INDEPENDENT-STRUCT-(1,1)'};rows={};
    for i=1:size(pairs,1)
        x=find(names==pairs{i,2});y=find(names==pairs{i,3});
        for m={'DeltaM','OmittedOffline','OmittedOracle','RE','Overlap'}
            for k=1:nK
                d=M.(m{1})(:,y,k)-M.(m{1})(:,x,k);bs=W*d;base=w0*M.(m{1})(:,x,k);
                rows(end+1,:)={string(pairs{i,1}),string(pairs{i,2}),string(pairs{i,3}),string(m{1}),rho(k),w0*d,prctile(bs,2.5),prctile(bs,97.5),(w0*d)/base}; %#ok<AGROW>
            end
        end
    end
    P=cell2table(rows,'VariableNames',{'Pair','First','Second','Metric','Rho','MeanDiff','Lo','Hi','RelativeToFirst'});writetable(P,fullfile(out,'r29d_causal_pairs.csv'));
    result=struct('classification',cls,'BER_gate_A',berA,'BER_gate_B',berB,'BER_GATE_PASS',berA||berB);
    r22a_io('json',fullfile(out,'r29d_gates_classification.json'),result);
end

function [c,v]=best_control(m,M,s,comp,k,w0)
    if m=="EtaPooled",val=arrayfun(@(r)(w0*M.CapturedOffline(:,r,k))/(w0*M.CapturedOracle(:,r,k)),comp);[v,b]=max(val);
    elseif m=="JaccardOverall",val=arrayfun(@(r)mean(s.jac(:,r,k)),comp);[v,b]=max(val);
    else,val=w0*M.(m)(:,comp,k);[v,b]=min(val);end
    c=comp(b);
end

function [est,lo,hi]=relative(m,M,s,t,c,k,W,w0)
    % Lower-is-better metrics: relative reduction (ctrl - t)/ctrl (ratio of means). Eta/Jaccard: absolute increase.
    switch m
        case 'EtaPooled'
            f=@(w)(w*M.CapturedOffline(:,t,k))./(w*M.CapturedOracle(:,t,k))-(w*M.CapturedOffline(:,c,k))./(w*M.CapturedOracle(:,c,k));
            est=f(w0);bs=f(W);
        case 'JaccardOverall'
            d=double(s.jac(:,t,k)-s.jac(:,c,k));est=mean(d);rng(292300000+k,'twister');bs=mean(d(randi(numel(d),numel(d),2000)),1)';   % pair bootstrap
        otherwise
            x=M.(m)(:,t,k);y=M.(m)(:,c,k);f=@(w)(w*y-w*x)./(w*y);est=f(w0);bs=f(W);
    end
    lo=prctile(bs,2.5);hi=prctile(bs,97.5);
end

function ok=velocity_sign(M,comp,t,k,w0)
    % R_E reduction vs the best comparator has the same (favourable) sign in >= 4 of 5 velocities.
    [c,~]=best_control('RE',M,[],comp,k,w0);n=size(M.RE,1)/5;d=reshape(M.RE(:,c,k)-M.RE(:,t,k),n,5);ok=nnz(mean(d,1)>0)>=4;
end

function n=runs(mask)
    % Longest run of adjacent true entries.
    n=0;c=0;for x=mask(:)',if x,c=c+1;n=max(n,c);else,c=0;end,end
end
