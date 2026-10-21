function varargout=r29_families(action,varargin)
% R2.9 representation families on the 24 (delay, M_tau) x 12 (Doppler, N_nu) finite grid.
% U = kron(U_nu, U_tau) (column-major vec, delay fast), U_axis = synthesis matrix (columns = basis functions).
% All design families satisfy W1 (current application constraint class); W2 is not imposed.
% Axis kinds (reduced W1 coordinates):
%   DY  (M_FB=2, basic depth L): x_j in S^1 via angles (L params); P_L fixed by W1 (SO(1) stabilizer).
%       Applied through the current periodic Mallat tree (J=2 on N_nu=12, J=3 on M_tau=24).
%   M4B (M_FB=4, basic depth L, reduced lossless): x_j in S^3 via hyperspherical angles (3 each);
%       P_L = diag(1,R3)*A(I) with R3 in SO(3) (3 params). Total 3L+3. Single M_FB=4 level.
%   M4S (M_FB=4, structured K_FB, gauge-normal): W_0..W_{K-1} in SO(4) (6 each), phi_k, psi_k,
%       terminal W_K = diag(1,R3)*A(T) (3 params). Total 8K+3. Single M_FB=4 level.
%   Tree 'HYBRID' (delay axis only): single M_FB=4 level, then the frozen current delay-axis
%       dyadic bank (alpha=pi/4) on the 6-sample lowpass block (6 -> 3).
    switch action
        case 'spec',varargout{1}=spec(varargin{:});
        case 'nparams',varargout{1}=nparams(varargin{:});
        case 'axis_E',varargout{1}=axis_E(varargin{:});
        case 'axis_U',varargout{1}=axis_U(varargin{:});
        case 'U',[varargout{1:2}]=build_U(varargin{:});
        case 'transform',varargout{1}=transform(varargin{:});
        case 'objective',[varargout{1:2}]=objective(varargin{:});
        case 'tree',varargout{1}=mallat_tree(varargin{:});
        otherwise,error('R29F:Action','Unknown action %s.',action);
    end
end

function f=spec(name)
    ax=@(kind,depth,tree)struct('kind',kind,'depth',depth,'tree',tree);
    switch name
        case 'DYADIC-W1-D1',f=struct('name',name,'axes',{{ax('DY',1,'MALLAT'),ax('DY',1,'MALLAT')}},'shared',false);
        case 'DYADIC-W1-D3',f=struct('name',name,'axes',{{ax('DY',3,'MALLAT'),ax('DY',3,'MALLAT')}},'shared',false);
        case 'DYADIC-W1-D7',f=struct('name',name,'axes',{{ax('DY',7,'MALLAT'),ax('DY',7,'MALLAT')}},'shared',false);
        case 'M4-BASIC-D1',f=struct('name',name,'axes',{{ax('M4B',1,'SINGLE'),ax('M4B',1,'SINGLE')}},'shared',false);
        case 'M4-INDEPENDENT-STRUCT-(1,1)',f=struct('name',name,'axes',{{ax('M4S',1,'SINGLE'),ax('M4S',1,'SINGLE')}},'shared',false);
        case 'M4-SHARED-STRUCT-(1,1)',f=struct('name',name,'axes',{{ax('M4S',1,'SINGLE'),ax('M4S',1,'SINGLE')}},'shared',true);
        case 'M4-BASIC-(3,3)',f=struct('name',name,'axes',{{ax('M4B',3,'SINGLE'),ax('M4B',3,'SINGLE')}},'shared',false);
        case 'M4-INDEPENDENT-STRUCT-(1,2)',f=struct('name',name,'axes',{{ax('M4S',1,'SINGLE'),ax('M4S',2,'SINGLE')}},'shared',false);
        case 'M4-BASIC-(3,6)',f=struct('name',name,'axes',{{ax('M4B',3,'SINGLE'),ax('M4B',6,'SINGLE')}},'shared',false);
        case 'HYBRID-M4-DYADIC',f=struct('name',name,'axes',{{ax('M4S',1,'SINGLE'),ax('M4S',1,'HYBRID')}},'shared',false);
        otherwise,error('R29F:Family','Unknown family %s.',name);
    end
end

function n=axis_nparams(a)
    switch a.kind,case 'DY',n=a.depth;case 'M4B',n=3*a.depth+3;case 'M4S',n=8*a.depth+3;end
end

function n=nparams(f)
    if f.shared,n=axis_nparams(f.axes{1});else,n=axis_nparams(f.axes{1})+axis_nparams(f.axes{2});end
end

function x=hyper(c)
    % S^3 hyperspherical coordinates.
    x=[cos(c(1));sin(c(1))*cos(c(2));sin(c(1))*sin(c(2))*cos(c(3));sin(c(1))*sin(c(2))*sin(c(3))];
end

function E=axis_E(a,c)
    switch a.kind
        case 'DY'
            L=a.depth;x=[cos(c(:).');sin(c(:).')];P=r29_pu('w1_terminal',eye(1),eye(2));   % SO(1) stabilizer is trivial
            E=r29_pu('lossless',x,P);
        case 'M4B'
            L=a.depth;x=zeros(4,L);for j=1:L,x(:,j)=hyper(c(3*j-2:3*j));end
            P=r29_pu('w1_terminal',r29_pu('so',3,c(3*L+1:3*L+3)),eye(4));E=r29_pu('lossless',x,P);
        case 'M4S'
            K=a.depth;W=cell(1,K+1);for k=1:K,W{k}=r29_pu('so',4,c(6*(k-1)+(1:6)));end;W{K+1}=eye(4);
            phi=c(6*K+(1:K));psi=c(7*K+(1:K));T=r29_pu('E1',r29_pu('struct_normal',W,phi,psi));
            W{K+1}=r29_pu('w1_terminal',r29_pu('so',3,c(8*K+1:8*K+3)),T);E=r29_pu('struct_normal',W,phi,psi);
    end
end

function U=mallat_tree(E,N)
    if mod(N,2)==1,U=eye(N);return,end
    S=r29_pu('block_matrix',E,N)';half=N/2;U=S*blkdiag(mallat_tree(E,half),eye(half));
end

function U=axis_U(a,c,N)
    persistent P6
    if isempty(P6),[P6,~,~]=fir_mallat_r22a(6,pi/4);end   % frozen current delay-axis dyadic bank (alpha_tau = pi/4)
    E=axis_E(a,c);
    switch a.tree
        case 'MALLAT',U=mallat_tree(E,N);
        case 'SINGLE',U=r29_pu('block_matrix',E,N)';
        case 'HYBRID'
            B=N/4;U=r29_pu('block_matrix',E,N)'*blkdiag(P6,eye(N-B));assert(B==6);
    end
end

function [U,axes]=build_U(f,c)
    c=c(:).';n1=axis_nparams(f.axes{1});
    if f.shared,c1=c;c2=c;else,c1=c(1:n1);c2=c(n1+1:end);end
    Unu=axis_U(f.axes{1},c1,12);Utau=axis_U(f.axes{2},c2,24);U=kron(Unu,Utau);axes={Unu,Utau};
end

function H=transform(G,Unu,Utau)
    % H = U.' G U for U = kron(Unu,Utau) (real), via mode products; G: 288 x 288 x C.
    C=size(G,3);X=reshape(G,24,12,[]);X=pagemtimes(pagemtimes(Utau.',X),Unu);X=reshape(X,288,288,C);
    X=pagetranspose(X);Y=reshape(X,24,12,[]);Y=pagemtimes(pagemtimes(Utau.',Y),Unu);H=pagetranspose(reshape(Y,288,288,C));
end

function [J,e]=objective(f,c,G,rhos)
    % e(U,rho) = mean over G of ||H - M.*H||^2/||H||^2 with M = top-K of mean |H|^2 over the same bank.
    [~,ax]=build_U(f,c);C=size(G,3);sup=zeros(288);nrm=zeros(288);chunk=50;
    for s=1:chunk:C
        idx=s:min(C,s+chunk-1);a=abs(transform(G(:,:,idx),ax{1},ax{2})).^2;tot=sum(a,[1,2]);
        sup=sup+sum(a,3);nrm=nrm+sum(a./tot,3);
    end
    sup=sup/C;nrm=nrm/C;[~,order]=sort(sup(:),'descend');cum=cumsum(nrm(order));
    e=max(1-cum(round(rhos*288^2)).',1e-16);J=mean(log10(e));
end
