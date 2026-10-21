function result=r29_stage_refine(paths,p,mode,families)
% Protocol revision R2.9c-rev1 (DESIGN only), applied identically to every family.
% Validated issue: the top-K objective is discontinuous (support switching); fminunc stops at the first step
% (exit flag 5) and random screening misses isolated special points: DYADIC-W1-D1 reached J = -2.637 while its
% own 2-D grid minimum (J = -2.8721) is the current Independent PUFB point (theta_nu = pi, theta_tau = pi/4).
% Revision (the preregistered fallback "pattern search", implemented locally; Global Optimization Toolbox absent):
%   screening  : previous 64 random points + 64 lattice points (every coordinate a multiple of pi/4)
%   starts     : 3 best distinct previous finals (pattern search only)
%              + 5 best new screening points not previously used (fminunc as before, then pattern search)
%   polish     : generalized pattern search, opportunistic coordinate polling, initial step pi/8, halving to 1e-4,
%                <= 100*(d+2) evaluations per start
% Previous runs are kept; VALIDATION selection uses the union of both runs.
    if nargin<3,mode='BROAD';end
    if nargin<4,families=p.families_all(~strcmp(p.families_all,'HYBRID-M4-DYADIC'));end
    out=fullfile(paths.c,'data');G=r29_io('load_split',paths,'design');
    Gc=parallel.pool.Constant(@()r29_io('load_split',paths,'design'));
    switch mode
        case 'BROAD'
            for i=1:numel(families)
                name=families{i};file=fullfile(out,['r29c_refine_',r29_io('tag',name),'.mat']);if isfile(file),continue,end
                old=load(fullfile(out,['r29c_broad_',r29_io('tag',name),'.mat']));
                [runs,screen]=refine_family(name,G,Gc,p.rho_design,old,sum(double(name))+293400); %#ok<ASGLU>
                save(file,'runs','screen','name');fprintf('R29 REFINE %s best J=%.5f (previous %.5f)\n',name,min(runs.JFinal),min(old.runs.JFinal));
            end
        case 'BUDGET'
            for name=["M4-INDEPENDENT-STRUCT-(1,1)","M4-BASIC-(3,3)"]
                for rho=p.rho_diag
                    file=fullfile(out,sprintf('r29c_refine_budget_%s_rho%03d.mat',r29_io('tag',char(name)),round(100*rho)));if isfile(file),continue,end
                    old=load(fullfile(out,sprintf('r29c_budget_%s_rho%03d.mat',r29_io('tag',char(name)),round(100*rho))));
                    ref=load(fullfile(out,['r29c_refine_',r29_io('tag',char(name)),'.mat']));   % refined broad finals as extra starts
                    [~,o]=sort(ref.runs.JFinal);extra=ref.runs.Params(o(1:2),:);
                    runs=polish_only(char(name),Gc,rho,[old.runs.Params;extra],["FROM_BUDGET";"FROM_BUDGET";"FROM_BUDGET";"FROM_REFINED_BROAD";"FROM_REFINED_BROAD"]); %#ok<NASGU>
                    save(file,'runs','name','rho');fprintf('R29 REFINE BUDGET %s rho=%.2f best=%.5f\n',name,rho,min(runs.JFinal));
                end
            end
    end
    result=true;
end

function [runs,screen]=refine_family(name,G,Gc,rhos,old,seed)
    f=r29_families('spec',name);d=r29_families('nparams',f);obj=@(c)r29_families('objective',f,c,G,rhos);rng(seed,'twister');
    L=(randi(8,64,d)-4)*pi/4;Jl=zeros(64,1);for i=1:64,Jl(i)=obj(L(i,:));end
    screen=table([old.screen.J;Jl],[repmat("RANDOM",64,1);repmat("LATTICE",64,1)],'VariableNames',{'J','Kind'});screen.Params=[old.screen.Params;L];
    used=old.runs.Params;[~,o]=sort(screen.J);fresh=zeros(0,d);
    for i=o',x=screen.Params(i,:);if min(vecnorm(used-x,2,2))>1e-9&&size(fresh,1)<5,fresh(end+1,:)=x;end,end %#ok<AGROW>
    [~,o]=sort(old.runs.JFinal);prev=old.runs.Params(o(1:3),:);
    starts=[prev;fresh];kind=[repmat("PREVIOUS_FINAL",3,1);repmat("NEW_SCREEN",size(fresh,1),1)];
    runs=run_starts(f,Gc,rhos,starts,kind);
end

function runs=polish_only(name,Gc,rhos,starts,kind)
    runs=run_starts(r29_families('spec',name),Gc,rhos,starts,kind(1:size(starts,1)));runs.Init(:)=kind(1:size(starts,1));
end

function runs=run_starts(f,Gc,rhos,starts,kind)
    n=size(starts,1);d=size(starts,2);Jinit=zeros(n,1);Jfin=Jinit;evals=Jinit;secs=Jinit;P=zeros(n,d);E=zeros(n,numel(rhos));
    opts=optimoptions('fminunc','Algorithm','quasi-newton','FiniteDifferenceType','forward','FiniteDifferenceStepSize',1e-6, ...
        'MaxIterations',150,'MaxFunctionEvaluations',150*(d+2),'FunctionTolerance',1e-6,'OptimalityTolerance',1e-8,'StepTolerance',1e-10,'Display','off');
    parfor s=1:n
        o=@(c)r29_families('objective',f,c,Gc.Value,rhos);t=tic;c=starts(s,:);J=o(c);Jinit(s)=J;ev=1;
        if kind(s)=="NEW_SCREEN",try,[c2,J2,~,out]=fminunc(o,c,opts);ev=ev+out.funcCount;if J2<J,c=c2;J=J2;end,catch,end,end
        [c,J,e2]=gps(o,c,J,pi/8,1e-4,100*(d+2));ev=ev+e2;
        [~,e]=o(c);Jfin(s)=J;P(s,:)=c;evals(s)=ev;secs(s)=toc(t);E(s,:)=e;
    end
    runs=table((1:n)',kind,Jinit,Jfin,evals,secs,'VariableNames',{'Start','Init','JInit','JFinal','Evaluations','Seconds'});runs.Params=P;runs.eDesign=E;
end

function [x,f,ev]=gps(obj,x,f,step,minstep,maxev)
    % Generalized pattern search on the coordinate basis (opportunistic polling, step halving on an unsuccessful poll).
    ev=0;d=numel(x);
    while step>=minstep&&ev<maxev
        improved=false;
        for i=1:d
            for sgn=[1,-1]
                y=x;y(i)=y(i)+sgn*step;fy=obj(y);ev=ev+1;
                if fy<f-1e-12,x=y;f=fy;improved=true;break,end
            end
            if ev>=maxev,break,end
        end
        if ~improved,step=step/2;end
    end
end
