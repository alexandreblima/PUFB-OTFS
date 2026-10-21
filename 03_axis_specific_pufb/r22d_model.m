function result=r22d_model(action,varargin)
    switch action
        case 'transform'
            result=transform(varargin{:});
        case 'metric'
            [sparse_mse,full_mse]=metric(varargin{:});result={sparse_mse,full_mse};
        case 'score'
            result=score(varargin{:});
        case 'support'
            result=support(varargin{:});
        case 'read'
            result=readtable(varargin{1},'Delimiter',',','VariableNamingRule','preserve','TextType','string');
        case 'hashes'
            paths=varargin{1};hashes=strings(numel(paths),1);
            for index=1:numel(paths),hashes(index)=string(r22a_io('hash',paths{index}));end
            result=table(string(paths(:)),hashes,'VariableNames',{'AbsolutePath','SHA256'});
        case 'selftest'
            result=selftest();
        otherwise
            error('R22D:Action','Unknown model operation.');
    end
end

function method=transform(alpha_nu,alpha_tau)
    assert(all([alpha_nu,alpha_tau]>=0&[alpha_nu,alpha_tau]<pi),'R22D:Domain','Canonical domain [0,pi).');
    [Pnu,~,nu]=fir_mallat_r22a(12,alpha_nu);[Ptau,~,tau]=fir_mallat_r22a(24,alpha_tau);U=kron(Pnu,Ptau);
    residuals=[norm(Pnu'*Pnu-eye(12),'fro'),norm(Ptau'*Ptau-eye(24),'fro'),norm(U'*U-eye(288),'fro'), ...
        norm(U*U'-eye(288),'fro'),nu.ParaunitarityResidual,tau.ParaunitarityResidual, ...
        nu.DCResidual,tau.DCResidual,nu.HighpassDCResidual,tau.HighpassDCResidual];
    assert(max(residuals)<1e-11,'R22D:StructuralGate','PR/paraunitarity/DC/vanishing moment failed.');
    assert(isreal(U)&&nu.FilterLength==4&&tau.FilterLength==4&&nu.Depth==2&&tau.Depth==3 ...
        &&isequal(size(nu.E),[2,2,2])&&isequal(size(tau.E),[2,2,2]),'R22D:Architecture','Architecture changed.');
    method=struct('Name',"",'AlphaNu',alpha_nu,'AlphaTau',alpha_tau,'Pnu',Pnu,'Ptau',Ptau, ...
        'U',U,'Nu',nu,'Tau',tau,'Residuals',residuals,'Mask',[]);
end

function [sparse_mse,full_mse]=metric(H,mask,snr_db,full_singular)
    Hs=H.*mask;gram=Hs'*Hs;
    [vectors,eigenvalues]=eig((gram+gram')/2,'vector');eigenvalues=max(real(eigenvalues),0);
    projected=vectors'*Hs'*H;row_energy=sum(abs(projected).^2,2);coupling=real(diag(projected*vectors));
    dimension=size(H,1);sparse_mse=zeros(1,numel(snr_db));full_mse=sparse_mse;
    for point=1:numel(snr_db)
        variance=10^(-snr_db(point)/10);denominator=eigenvalues+variance;
        sparse_mse(point)=(dimension-2*sum(coupling./denominator)+sum(row_energy./denominator.^2) ...
            +variance*sum(eigenvalues./denominator.^2))/dimension;
        full_mse(point)=mean(variance./(full_singular.^2+variance));
    end
    assert(all(isfinite(sparse_mse))&&all(sparse_mse>=full_mse-1e-10),'R22D:MSE','Invalid exact mismatch MSE.');
end

function mask=support(bank,U,K)
    energy=zeros(288);
    for channel=1:size(bank,3),H=U'*bank(:,:,channel)*U;energy=energy+abs(H).^2;end
    energy=energy/size(bank,3);[~,order]=sort(energy(:),'descend');mask=false(288);mask(order(1:K))=true;
    assert(nnz(mask)==K,'R22D:K','Support cardinality changed.');
end

function result=score(pair,bank,singular,baseline,snr_db,p,frozen_mask)
    method=transform(pair(1),pair(2));
    if isempty(frozen_mask),method.Mask=support(bank,method.U,p.K);else,method.Mask=frozen_mask;end
    assert(nnz(method.Mask)==p.K);delta=zeros(size(bank,3),numel(snr_db));
    for channel=1:size(bank,3)
        H=method.U'*bank(:,:,channel)*method.U;
        [sparse_mse,full_mse]=metric(H,method.Mask,snr_db,singular(:,channel));delta(channel,:)=sparse_mse-full_mse;
    end
    strata={p.strata.LOW,p.strata.MID,p.strata.HIGH};raw=zeros(1,3);
    for stratum=1:3,raw(stratum)=mean(delta(:,ismember(snr_db,strata{stratum})),'all');end
    if isempty(baseline),ratios=ones(1,3);else,ratios=raw./baseline;end
    result=struct('Pair',pair,'Jrx',mean(ratios),'Ratios',ratios,'Raw',raw,'Delta',delta,'Method',method);
end

function result=selftest()
    pairs=[pi/3,pi/3;.78797323934301,.7967306572;.1,2.9;0,pi/2];worst=0;
    stream=RandStream('mt19937ar','Seed',119001);H=randn(stream,12)+1i*randn(stream,12);
    masks={abs(H)>1,true(12),false(12)};
    for index=1:size(pairs,1),method=transform(pairs(index,1),pairs(index,2));assert(max(method.Residuals)<1e-11);end
    for index=1:numel(masks)
        [calculated,full]=metric(H,masks{index},3:3:24,svd(H));
        [direct,direct_full]=confirm_exact_mse(H,masks{index},3:3:24);
        worst=max([worst,max(abs(calculated-direct)),max(abs(full-direct_full))]);
    end
    assert(worst<1e-11,'R22D:MetricCheck','Spectral solve disagrees with direct Cholesky residual.');
    result=struct('PASS',true,'ExactMetricDirectSolveError',worst,'PhysicalChannelsGenerated',false);
end