function result=r29_stage_design(paths,p,mode,families)
% DESIGN-only optimization. Objective J_broad = mean_rho log10 e(U,rho), rho in {0.3,...,0.7},
% DESIGN mean-energy support (approved objective; the bounded ratio score is a VALIDATION diagnostic only).
% Starts per family: 64-point random screening; 7 best screened points + 1 structured (Haar/D4-like) seed.
% Optimizer: fminunc quasi-Newton, forward FD gradient, <=150 iterations, relative function tolerance 1e-6.
% Fallback (documented deviation): fminsearch, because patternsearch (Global Optimization Toolbox) is not installed.
    if nargin<3,mode='BROAD';end
    if nargin<4,families=p.families;end
    G=load_split(paths,'design');out=fullfile(paths.c,'data');
    if ~strcmp(mode,'TIME'),Gc=parallel.pool.Constant(@()r29_io('load_split',paths,'design'));end   % workers read DESIGN from disk (no 400 MB broadcast)
    switch mode
        case 'TIME'
            f=r29_families('spec','M4-INDEPENDENT-STRUCT-(1,1)');c=rand(1,r29_families('nparams',f));
            tic;for i=1:3,r29_families('objective',f,c,G,p.rho_design);end;result=struct('seconds_per_evaluation',toc/3);disp(result);return
        case 'BROAD'
            for i=1:numel(families)
                name=families{i};file=fullfile(out,['r29c_broad_',tag(name),'.mat']);if isfile(file),continue,end
                [runs,screen]=optimize_family(name,G,Gc,p.rho_design,[],sum(double(name))+293100); %#ok<ASGLU>
                save(file,'runs','screen','name');fprintf('R29 DESIGN %s best J=%.5f\n',name,min(runs.JFinal));
            end
        case 'BUDGET'
            for name=["M4-INDEPENDENT-STRUCT-(1,1)","M4-BASIC-(3,3)"]
                broad=load(fullfile(out,['r29c_broad_',tag(char(name)),'.mat']));
                for rho=p.rho_diag
                    file=fullfile(out,sprintf('r29c_budget_%s_rho%03d.mat',tag(char(name)),round(100*rho)));if isfile(file),continue,end
                    [~,o]=sort(broad.runs.JFinal);starts=broad.runs.Params(o(1:3),:);   % broad optimum + next two distinct finals
                    [runs,screen]=optimize_family(char(name),G,Gc,rho,starts,round(1000*rho)+293200); %#ok<ASGLU>
                    save(file,'runs','screen','name','rho');fprintf('R29 BUDGET %s rho=%.2f best=%.5f\n',name,rho,min(runs.JFinal));
                end
            end
    end
    result=summarize(out);
end

function G=load_split(paths,split)
    T=readtable(fullfile(paths.data,sprintf('r29_%s_manifest.csv',split)),'TextType','string','Delimiter',',','VariableNamingRule','preserve');
    seal=jsondecode(fileread(fullfile(paths.data,sprintf('r29_%s_seal.json',split))));
    assert(strcmp(seal.manifest_sha256,r22a_io('hash',fullfile(paths.data,sprintf('r29_%s_manifest.csv',split)))),'manifest changed');
    G=complex(zeros(288,288,height(T)));
    for i=1:height(T),l=load(T.AbsolutePath(i),'G');G(:,:,i)=l.G;end
end

function s=tag(name)
    s=regexprep(name,'[^A-Za-z0-9]+','_');
end

function [runs,screen]=optimize_family(name,G,Gc,rhos,starts,seed)
    f=r29_families('spec',name);d=r29_families('nparams',f);obj=@(c)r29_families('objective',f,c,G,rhos);
    rng(seed,'twister');
    if isempty(starts)
        X=random_points(f,64);Js=zeros(64,1);for i=1:64,Js(i)=obj(X(i,:));end
        [~,o]=sort(Js);starts=[X(o(1:7),:);structured_seed(f)];kind=[repmat("SCREEN",7,1);"STRUCTURED"];
        screen=table(Js,'VariableNames',{'J'});screen.Params=X;
    else
        kind=repmat("FROM_BROAD",size(starts,1),1);screen=table();
    end
    n=size(starts,1);Jinit=zeros(n,1);Jfin=zeros(n,1);evals=zeros(n,1);iters=zeros(n,1);flag=zeros(n,1);fallback=false(n,1);P=zeros(n,d);secs=zeros(n,1);
    opts=optimoptions('fminunc','Algorithm','quasi-newton','FiniteDifferenceType','forward','FiniteDifferenceStepSize',1e-6, ...
        'MaxIterations',150,'MaxFunctionEvaluations',150*(d+2),'FunctionTolerance',1e-6,'OptimalityTolerance',1e-8,'StepTolerance',1e-10,'Display','off');
    parfor s=1:n
        objw=@(c)r29_families('objective',f,c,Gc.Value,rhos);t=tic;c0=starts(s,:);Jinit(s)=objw(c0);
        try
            [c,J,ef,o]=fminunc(objw,c0,opts);ev=o.funcCount;it=o.iterations;fb=false;
        catch
            ef=-99;J=Jinit(s);c=c0;ev=0;it=0;fb=true;
        end
        if ef<=0||J>=Jinit(s)-1e-9   % no admissible progress -> derivative-free fallback from the same start
            [c2,J2,ef2,o2]=fminsearch(objw,c0,optimset('MaxIter',150*d,'MaxFunEvals',150*(d+2),'TolFun',1e-6,'TolX',1e-8,'Display','off'));
            ev=ev+o2.funcCount;it=it+o2.iterations;fb=true;if J2<J,c=c2;J=J2;ef=ef2;end
        end
        Jfin(s)=J;P(s,:)=c;evals(s)=ev;iters(s)=it;flag(s)=ef;fallback(s)=fb;secs(s)=toc(t);
    end
    [~,e]=cellfun(@(c)obj(c),num2cell(P,2),'UniformOutput',false);E=cell2mat(e);
    runs=table((1:n)',kind,Jinit,Jfin,evals,iters,flag,fallback,secs,'VariableNames',{'Start','Init','JInit','JFinal','Evaluations','Iterations','ExitFlag','Fallback','Seconds'});
    runs.Params=P;runs.eDesign=E;
end

function X=random_points(f,n)
    % Uniform angles; polar hyperspherical angles drawn away from the poles; phi,psi in (0,pi).
    d=r29_families('nparams',f);X=2*pi*rand(n,d)-pi;
    for a=axis_slices(f)'
        ax=f.axes{a.axis};idx=a.offset+(1:a.n);
        switch ax.kind
            case 'M4B',for j=1:ax.depth,X(:,idx(3*j-2:3*j-1))=0.2+(pi-0.4)*rand(n,2);end
            case 'M4S',K=ax.depth;X(:,idx(6*K+(1:2*K)))=0.1+(pi-0.2)*rand(n,2*K);
                X(:,idx(1:6*K))=pi*(rand(n,6*K)-0.5);X(:,idx(8*K+(1:3)))=pi*(rand(n,3)-0.5);
        end
    end
end

function s=axis_slices(f)
    n1=nparams_axis(f.axes{1});
    if f.shared,s=struct('axis',{1},'offset',{0},'n',{n1});
    else,s=struct('axis',{1,2},'offset',{0,n1},'n',{n1,nparams_axis(f.axes{2})});end
    s=s(:);
end

function n=nparams_axis(a)
    switch a.kind,case 'DY',n=a.depth;case 'M4B',n=3*a.depth+3;case 'M4S',n=8*a.depth+3;end
end

function c=structured_seed(f)
    % DY: D4 first factor (current alpha = pi/3 locus) followed by pure delays (x = e1).
    % M4B/M4S: Walsh-Hadamard-like point (delays only / identity rotations); small offsets avoid chart poles.
    c=zeros(1,r29_families('nparams',f));
    for a=axis_slices(f)'
        ax=f.axes{a.axis};idx=a.offset+(1:a.n);
        switch ax.kind
            case 'DY'
                [~,~,meta]=fir_mallat_r22a(6,pi/3);E1=[meta.h(3),meta.h(4);meta.g(3),meta.g(4)];x0=null(E1);c(idx(1))=atan2(x0(2),x0(1));
            case 'M4B',for j=1:ax.depth,c(idx(3*j-2:3*j))=[1e-2,pi/2,pi/2];end
            case 'M4S',K=ax.depth;c(idx(6*K+(1:2*K)))=1e-2;
        end
    end
end

function S=summarize(out)
    files=dir(fullfile(out,'r29c_*.mat'));rows={};
    for i=1:numel(files)
        l=load(fullfile(out,files(i).name));r=l.runs;
        rows(end+1,:)={string(files(i).name),height(r),min(r.JFinal),median(r.JFinal),sum(r.Evaluations),sum(r.Fallback),sum(r.ExitFlag<=0),sum(r.Seconds)}; %#ok<AGROW>
    end
    S=cell2table(rows,'VariableNames',{'File','Starts','BestJ','MedianJ','Evaluations','Fallbacks','NonpositiveExit','Seconds'});
    writetable(S,fullfile(out,'r29c_design_summary.csv'));disp(S);
end
