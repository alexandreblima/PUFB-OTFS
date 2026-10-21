function result=r29_stage_test(paths,p,mode)
% Frozen-TEST compression/support analysis (sec. 18-20). Requires the candidate-freeze manifest.
% ANALYTIC: R2.6 receiver decomposition (arm R, SNR 3:3:24 dB, lambda=sigma^2, no eta), offline (DESIGN) and oracle supports.
% SUPPORT : p_j maps, H_support, offline/oracle overlap, eta_capture, R_E, Jaccard (predeclared pairs), regime transfer.
    out=fullfile(paths.d,'data');freeze=fullfile(paths.c,'config','r29_candidate_freeze.mat');
    man=jsondecode(fileread(fullfile(paths.c,'config','r29_candidate_freeze_manifest.json')));
    assert(strcmp(man.freeze_mat_sha256,r22a_io('hash',freeze)),'freeze file changed after manifest');
    R=representations(paths,freeze);
    switch mode
        case 'SELFTEST',result=selftest(paths,p,R);
        case 'OPEN'
            opening=fullfile(paths.data,'r29_test_opening.json');assert(~isfile(opening));
            r22a_io('json',opening,struct('opened_utc',char(datetime('now','TimeZone','UTC')),'freeze_manifest_sha256', ...
                r22a_io('hash',fullfile(paths.c,'config','r29_candidate_freeze_manifest.json')),'representations',{{R.name}}));result=true;
        case 'ORDERS',result=orders(paths,p,R,out);
        case 'ANALYTIC',result=analytic(paths,p,R,out);
        case 'SUPPORT',result=support(paths,p,R,out);
    end
end

function R=representations(paths,freeze)
    C=r29_io('controls',paths);l=load(freeze,'cand');R=struct('name',{},'group',{},'A',{},'B',{});
    for c=C,R(end+1)=struct('name',c.name,'group','CONTROL','A',c.U,'B',[]);end %#ok<AGROW>
    for c=l.cand
        if strcmp(c.label,'BROAD'),name=c.name;else,name=[c.name,'@',c.label];end
        if startsWith(c.name,'DYADIC'),g='M2_NEW';elseif startsWith(c.name,'HYBRID'),g='HYBRID';else,g='M4';end
        if ~strcmp(c.label,'BROAD'),g='DIAG';end
        R(end+1)=struct('name',name,'group',g,'A',c.Unu,'B',c.Utau); %#ok<AGROW>
    end
end

%% ------------------------------------------------------------------ DESIGN supports (frozen before TEST metrics)
function result=orders(paths,p,R,out)
    file=fullfile(out,'r29d_offline_orders.mat');assert(~isfile(file));
    Gd=r29_io('load_split',paths,'design');T=readtable(fullfile(paths.data,'r29_design_manifest.csv'),'Delimiter',',','VariableNamingRule','preserve');
    n2=288^2;order=zeros(n2,numel(R),'uint32');order_v=zeros(n2,numel(R),5,'uint32');
    for r=1:numel(R)
        e=r29_io('mean_energy',R(r).A,R(r).B,Gd);[~,o]=sort(e(:),'descend');order(:,r)=o;
        for v=1:5,e=r29_io('mean_energy',R(r).A,R(r).B,Gd(:,:,T.VelocityIndex==v));[~,o]=sort(e(:),'descend');order_v(:,r,v)=o;end
    end
    names={R.name};save(file,'order','order_v','names');result=struct('representations',numel(R),'sha256',r22a_io('hash',file));
end

%% ------------------------------------------------------------------ ANALYTIC
function result=analytic(paths,p,R,out)
    file=fullfile(out,'r29d_analytic.mat');assert(~isfile(file));
    o=load(fullfile(out,'r29d_offline_orders.mat'),'order');order=o.order;
    T=readtable(fullfile(paths.data,'r29_test_manifest.csv'),'TextType','string','Delimiter',',','VariableNamingRule','preserve');
    files=T.AbsolutePath;K=round(p.rho_test*288^2);snr=3:3:24;nR=numel(R);count=height(T);
    A={R.A};B={R.B};direct=ismember(1:count,[1:3,201:203,401:403,601:603,801:803]);
    terms=zeros(count,nR,numel(K),numel(snr),5,2);   % T_reg T_omit T_cross T_noise M_true ; policy OFFLINE/ORACLE
    scalar=zeros(count,nR,numel(K),4,2);                      % omitted (normalized), receiver norm ||W||^2/n at 12 dB, rank, captured (normalized)
    full=zeros(count,numel(snr));residual=zeros(count,2);
    if isempty(gcp('nocreate')),parpool('Processes',6);end
    parfor c=1:count
        l=load(files(c),'G');[terms(c,:,:,:,:,:),scalar(c,:,:,:,:),full(c,:),residual(c,:)]=channel_analytic(l.G,A,B,order,K,snr,direct(c)); %#ok<PFBNS>
    end
    names={R.name};groups={R.group};save(file,'terms','scalar','full','residual','names','groups','K','snr','-v7.3');
    rho1=numel(K);dm=squeeze(double(terms(:,:,rho1,:,5,1)))-reshape(full,count,1,[]);
    result=struct('max_abs_DeltaM_rho1',max(abs(dm),[],'all'),'max_omitted_rho1',max(scalar(:,:,rho1,1,:),[],'all'), ...
        'decomposition_vs_formula',max(residual(:,1)),'decomposition_vs_direct_cholesky',max(residual(:,2)));
    result.PASS=result.max_abs_DeltaM_rho1<1e-10&&result.max_omitted_rho1<1e-12&&result.decomposition_vs_formula<1e-10&&result.decomposition_vs_direct_cholesky<1e-9;
    r22a_io('json',fullfile(out,'r29d_analytic_gate.json'),result);
end

function [terms,scalar,full,residual]=channel_analytic(G,A,B,order,K,snr,direct)
    n=288;nR=numel(A);nK=numel(K);nS=numel(snr);variance=10.^(-snr/10);
    terms=zeros(nR,nK,nS,5,2);scalar=zeros(nR,nK,4,2);s_full=svd(G);full=mean(variance./(s_full.^2+variance),1);residual=[0,0];
    for r=1:nR
        H=r29_io('coeffs',A{r},B{r},G);energy=abs(H).^2;total=sum(energy,'all');[~,oracle]=sort(energy(:),'descend');
        pol={double(order(:,r)),oracle};
        for policy=1:2
            for k=1:nK
                mask=false(n);mask(pol{policy}(1:K(k)))=true;Hs=H.*mask;E=H-Hs;
                [L,S,Rr]=svd(Hs);s=diag(S);LE=L'*E;row_energy=sum(abs(LE).^2,2);cross_diag=real(diag(LE*Rr));
                LH=L'*H;coupling=real(diag(LH*Rr));row_full=sum(abs(LH).^2,2);
                for j=1:nS
                    g=s./(s.^2+variance(j));gs=g.*s;
                    t=[sum((gs-1).^2),sum(g.^2.*row_energy),2*sum((gs-1).*g.*cross_diag),variance(j)*sum(g.^2)]/n;
                    formula=(n-2*sum(g.*coupling)+sum(g.^2.*row_full)+variance(j)*sum(g.^2))/n;
                    terms(r,k,j,:,policy)=[t,sum(t)];residual(1)=max(residual(1),abs(sum(t)-formula));
                    if direct&&policy==1&&ismember(j,[1,nS])&&ismember(k,[1,8,nK])
                        f=chol(Hs'*Hs+variance(j)*eye(n));W=f\(f'\Hs');Am=W*Hs-eye(n);Bm=W*E;
                        d=[norm(Am,'fro')^2,norm(Bm,'fro')^2,2*real(sum(Am(:).*conj(Bm(:)))),variance(j)*norm(W,'fro')^2]/n;
                        residual(2)=max([residual(2),max(abs(d-t)),abs((norm(W*H-eye(n),'fro')^2+variance(j)*norm(W,'fro')^2)/n-sum(t))]);
                    end
                end
                g12=s./(s.^2+10^(-1.2));scalar(r,k,:,policy)=[norm(E,'fro')^2/total,sum(g12.^2)/n,nnz(s>n*eps(max(s))),sum(energy(mask))/total];
            end
        end
    end
end

%% ------------------------------------------------------------------ SUPPORT
function result=support(paths,p,R,out)
    file=fullfile(out,'r29d_support.mat');assert(~isfile(file));
    o=load(fullfile(out,'r29d_offline_orders.mat'));
    T=readtable(fullfile(paths.data,'r29_test_manifest.csv'),'TextType','string','Delimiter',',','VariableNamingRule','preserve');
    count=height(T);vel=T.VelocityIndex;n2=288^2;K=round(p.rho_test*n2);nK=numel(K);nR=numel(R);
    % Predeclared pair samples (seed 292100000 overall; 292100000+v within velocity v).
    rng(292100000,'twister');pairs=zeros(2000,2);i=0;while i<2000,a=randi(count,1,2);if a(1)~=a(2),i=i+1;pairs(i,:)=a;end,end
    vp=cell(1,5);for v=1:5,rng(292100000+v,'twister');idx=find(vel==v);q=zeros(400,2);i=0;
        while i<400,a=idx(randi(numel(idx),1,2));if a(1)~=a(2),i=i+1;q(i,:)=a;end,end;vp{v}=q;end
    G=complex(zeros(288,288,count));for c=1:count,l=load(T.AbsolutePath(c),'G');G(:,:,c)=l.G;end
    p_map=zeros(n2,nR,nK,'single');Hs=zeros(nR,nK);eta=zeros(count,nR,nK);overlap=eta;RE=eta;jac=zeros(2000,nR,nK,'single');jacv=zeros(400,5,nR,nK,'single');
    transfer=zeros(5,5,nR,nK);   % rows: DESIGN velocity used for support; cols: TEST velocity; pooled captured / oracle captured
    for r=1:nR
        rank=zeros(n2,count,'uint32');capt_or=zeros(count,nK);capt_off=zeros(count,nK);capt_tr=zeros(count,5,nK);
        for s=1:100:count
            idx=s:min(count,s+99);a=abs(r29_io('coeffs',R(r).A,R(r).B,G(:,:,idx))).^2;a=reshape(a,n2,[]);a=a./sum(a,1);
            [sorted,or]=sort(a,1,'descend');cs=cumsum(sorted,1);
            for j=1:numel(idx)
                rank(or(:,j),idx(j))=1:n2;capt_or(idx(j),:)=cs(K,j);ca=cumsum(a(double(o.order(:,r)),j));capt_off(idx(j),:)=ca(K);
                for v=1:5,cv=cumsum(a(double(o.order_v(:,r,v)),j));capt_tr(idx(j),v,:)=cv(K);end
            end
        end
        off_rank=zeros(n2,1,'uint32');off_rank(o.order(:,r))=1:n2;
        for k=1:nK
            m=rank<=K(k);p_map(:,r,k)=mean(m,2);q=double(p_map(:,r,k))/sum(p_map(:,r,k));q=q(q>0);Hs(r,k)=-sum(q.*log(q))/log(n2);
            overlap(:,r,k)=sum(m&(off_rank<=K(k)),1)'/K(k);eta(:,r,k)=capt_off(:,k)./capt_or(:,k);RE(:,r,k)=capt_or(:,k)-capt_off(:,k);
            inter=sum(m(:,pairs(:,1))&m(:,pairs(:,2)),1);jac(:,r,k)=inter./(2*K(k)-inter);
            for v=1:5,inter=sum(m(:,vp{v}(:,1))&m(:,vp{v}(:,2)),1);jacv(:,v,r,k)=inter./(2*K(k)-inter);
                for w=1:5,transfer(v,w,r,k)=mean(capt_tr(vel==w,v,k))/mean(capt_or(vel==w,k));end,end
        end
        fprintf('R29 SUPPORT %s done\n',R(r).name);
    end
    names={R.name};groups={R.group};
    save(file,'p_map','Hs','eta','overlap','RE','jac','jacv','transfer','pairs','vp','names','groups','K','-v7.3');
    result=struct('representations',nR,'max_eta',max(eta,[],'all'),'min_RE',min(RE,[],'all'));
end

%% ------------------------------------------------------------------ SELFTEST (no TEST data)
function result=selftest(paths,p,R)
    % (1) The analytic kernel reproduces stored R2.6 ORACLE terms for the six controls on R2.6 TEST channel 1 (support-independent policy).
    % (2) Separable coefficients equal the dense U'GU for every frozen candidate.
    r26=load(fullfile(paths.v2,'results','r2_6_budget_continuation','data','r26_analytic_R.mat'),'terms','K','snr');
    l=load(fullfile(paths.v2,'results','r2_2d_axis_specific','test_channels','channel_0001.mat'),'G');
    A={R(1:6).A};B={R(1:6).B};dummy=repmat(uint32((1:288^2)'),1,6);
    terms=channel_analytic(l.G,A,B,dummy,r26.K,r26.snr,false);
    err=max(abs(terms(:,:,:,:,2)-squeeze(r26.terms(1,:,:,:,:,2))),[],'all');
    sep=0;Gt=randn(288,288,2)+1i*randn(288,288,2);
    for r=7:numel(R),U=kron(R(r).A,R(r).B);H=r29_io('coeffs',R(r).A,R(r).B,Gt);sep=max(sep,max(abs(H(:,:,2)-U'*Gt(:,:,2)*U),[],'all'));
        sep=max(sep,norm(U'*U-eye(288),'fro'));end
    result=struct('r26_oracle_terms_max_abs_error',err,'separable_vs_dense_and_unitarity',sep,'PASS',err<1e-12&&sep<1e-10);
    r22a_io('json',fullfile(paths.d,'data','r29d_selftest.json'),result);
end
