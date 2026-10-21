function result=r22e_math(action,varargin)
    switch action
        case 'geometry'
            H=varargin{1};mask=varargin{2};Hs=H.*mask;[left,diagonal,right]=svd(Hs,'econ');s=diag(diagonal);E=H-Hs;
            projected=left'*H;omitted=left'*E;
            result=struct('Singular',s,'Left',left,'Right',right,'H',H,'Hs',Hs,'E',E, ...
                'Projected',projected,'Omitted',omitted,'RowEnergy',sum(abs(projected).^2,2), ...
                'Coupling',real(diag(projected*right)),'Q',norm(Hs,'fro')^2/size(Hs,1));
        case 'risk'
            result=risk(varargin{:});
        case 'receiver'
            Hs=varargin{1};variance=varargin{2};eta=varargin{3};loading=variance+eta*norm(Hs,'fro')^2/size(Hs,1);
            factor=chol(Hs'*Hs+loading*eye(size(Hs,2)));result=factor\(factor'\Hs');
        case 'select'
            values=varargin{1};etas=varargin{2};minimum=min(values);eligible=values<=minimum+1e-4*max(abs(minimum),realmin);
            result=find(eligible,1,'first');assert(issorted(etas));
        case 'bootstrap'
            values=varargin{1};count=varargin{2};seed=varargin{3};result=r22dcb_math('bootstrap',values,count,size(values,1)/count,seed,10000);
        case 'read'
            result=r22dcb_math('read',varargin{1});
        case 'selftest'
            H=reshape(sin(1:144),12,12)+1i*reshape(cos(3:146),12,12)+eye(12);mask=mod(reshape(1:144,12,12),3)~=0;
            geometry=r22e_math('geometry',H,mask);etas=[0,1e-6,.01,1];snrs=[20,30,36];values=risk(geometry,etas,snrs,false);maximum_error=0;
            for eta_index=1:numel(etas)
                for point=1:numel(snrs)
                    variance=10^(-snrs(point)/10);W=r22e_math('receiver',geometry.Hs,variance,etas(eta_index));
                    direct=(norm(W*H-eye(12),'fro')^2+variance*norm(W,'fro')^2)/12;
                    maximum_error=max(maximum_error,abs(values(eta_index,point)-direct));
                end
            end
            assert(maximum_error<1e-10&&r22e_math('select',[1,1.00001,2],[0,.1,1])==1);
            result=struct('PASS',true,'ExactRiskDirectSolveError',maximum_error,'QUsesSparseChannelOnly',true);
        otherwise
            error('R22E:Action','Unknown operation.');
    end
end

function result=risk(geometry,etas,snr_db,diagnostics)
    s=geometry.Singular;dimension=size(geometry.H,1);snr_db=snr_db(:).';etas=etas(:).';
    if diagnostics,result=zeros(numel(etas),numel(snr_db),11);else,result=zeros(numel(etas),numel(snr_db));end
    for eta_index=1:numel(etas)
        for point=1:numel(snr_db)
            variance=10^(-snr_db(point)/10);loading=variance+etas(eta_index)*geometry.Q;g=s./(s.^2+loading);
            mse=(dimension-2*sum(g.*geometry.Coupling)+sum(g.^2.*geometry.RowEnergy)+variance*sum(g.^2))/dimension;
            if ~diagnostics,result(eta_index,point)=mse;continue;end
            W=geometry.Right*(g.*geometry.Left');reg=W*geometry.Hs-eye(dimension);omit=W*geometry.E;
            reg_term=norm(reg,'fro')^2/dimension;omit_term=norm(omit,'fro')^2/dimension;
            cross=2*real(sum(reg(:).*conj(omit(:))))/dimension;noise=variance*norm(W,'fro')^2/dimension;
            assert(abs(mse-reg_term-omit_term-cross-noise)<1e-10*(1+mse));
            result(eta_index,point,:)=[mse,reg_term,omit_term,cross,noise,norm(W,'fro'),norm(omit,'fro'), ...
                (s(1)^2+loading)/(s(end)^2+loading),loading,geometry.Q,omit_term+cross];
        end
    end
end