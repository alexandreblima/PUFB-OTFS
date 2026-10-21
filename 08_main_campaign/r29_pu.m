function varargout=r29_pu(action,varargin)
% Paraunitary toolkit for R2.9, written from the definitions of
% geometry_constrained_learned_wavelets_rev16 (Sec. III-IV) and its supplement.
% Polyphase matrices are stored as M_FB x M_FB x (L+1) coefficient arrays: E(:,:,l+1) = coefficient of z^-l.
% Conventions:
%   Lambda(z) = diag(1, z^-1, ..., z^-1)                               (eq. intro-basic-cascade)
%   R_pq(phi) = expm(phi*J_pq), J_pq = e_p e_q' - e_q e_p', R_qp(phi) = R_pq(-phi)  (eqs. so-generators, givens-exp)
%   basic:      E = U_L Lambda U_{L-1} ... U_1 Lambda U_0
%   structured: L_k = W_k Lambda D'_k Lambda D_k Lambda, E = L_K...L_1 W_0, D_k = R_12(phi_k),
%               D'_k = G_k Q_k, G_k = R_{M,1}(psi_k) (ordered wrap plane), Q_k = internal shifted sections
%   Type-I row filters: H(z) = E(z^M) d(z), d(z) = (1, z^-1, ..., z^-(M-1))'  =>  h_m[M*l+p] = E_mp[l]
%   W1: P_L v = e_1 with v = 1/sqrt(M) (eq. compact-w1); here P_L = E(1).
    switch action
        case 'rot',varargout{1}=rot(varargin{:});
        case 'so',varargout{1}=so_from_coords(varargin{:});
        case 'basic',varargout{1}=basic(varargin{:});
        case 'lossless',varargout{1}=lossless(varargin{:});
        case 'struct_raw',varargout{1}=struct_raw(varargin{:});
        case 'struct_normal',[varargout{1:2}]=struct_normal(varargin{:});
        case 'partial_x',varargout{1}=partial_x(varargin{:});
        case 'gram',varargout{1}=gram(varargin{:});
        case 'filters',varargout{1}=filters(varargin{:});
        case 'E1',varargout{1}=sum(varargin{1},3);
        case 'paraunitarity',varargout{1}=paraunitarity(varargin{:});
        case 'block_matrix',varargout{1}=block_matrix(varargin{:});
        case 'apply_polyphase',varargout{1}=apply_polyphase(varargin{:});
        case 'apply_filterbank',varargout{1}=apply_filterbank(varargin{:});
        case 'polymul',varargout{1}=polymul(varargin{:});
        case 'w1_terminal',varargout{1}=w1_terminal(varargin{:});
        otherwise,error('R29:Action','Unknown action %s.',action);
    end
end

function R=rot(M,p,q,phi)
    J=zeros(M);J(p,q)=1;J(q,p)=-1;R=expm(phi*J);   % R_pq(phi); ordered pair (q,p) obtained by passing -phi
end

function U=so_from_coords(M,theta)
    % SO(M) chart: expm of the skew matrix with upper-triangular coordinates theta (length M(M-1)/2).
    A=zeros(M);A(triu(true(M),1))=theta;U=expm(A-A');
end

function Lam=delay(M)
    Lam=zeros(M,M,2);Lam(1,1,1)=1;for i=2:M,Lam(i,i,2)=1;end
end

function C=polymul(A,B)
    % (A*B)(z) for coefficient arrays.
    M=size(A,1);La=size(A,3);Lb=size(B,3);C=zeros(M,size(B,2),La+Lb-1);
    for i=1:La,for j=1:Lb,C(:,:,i+j-1)=C(:,:,i+j-1)+A(:,:,i)*B(:,:,j);end,end
end

function E=basic(U)
    % U: cell {U_0,...,U_L}. E = U_L Lambda ... U_1 Lambda U_0.
    M=size(U{1},1);E=reshape(U{1},M,M,1);Lam=delay(M);
    for k=2:numel(U),E=polymul(reshape(U{k},M,M,1),polymul(Lam,E));end
end

function E=lossless(x,P)
    % Reduced projector form (Prop. classical-lossless-cascade): E = P_L C_{x_{L-1}} ... C_{x_0},
    % C_x(z) = x x' + z^-1 (I - x x'). x: M x L matrix of unit vectors (columns x_0..x_{L-1}).
    M=size(P,1);E=reshape(eye(M),M,M,1);
    for j=1:size(x,2)
        v=x(:,j)/norm(x(:,j));C=zeros(M,M,2);C(:,:,1)=v*v';C(:,:,2)=eye(M)-v*v';E=polymul(C,E);
    end
    E=polymul(reshape(P,M,M,1),E);
end

function [E,U]=struct_normal(W,phi,psi)
    % Gauge-normal structured word (eq. gauge-normal-form):
    % E = W_K Lambda G_K Lambda D_K Lambda ... W_1 Lambda G_1 Lambda D_1 Lambda W_0.
    % W: cell {W_0,...,W_K}. Returns E and the equivalent basic factor list (depth 3K) for containment tests.
    M=size(W{1},1);K=numel(W)-1;U={W{1}};
    for k=1:K
        U{end+1}=rot(M,1,2,phi(k));U{end+1}=rot(M,1,M,-psi(k));U{end+1}=W{k+1}; %#ok<AGROW>  % R_{M,1}(psi)=R_{1,M}(-psi)
    end
    E=basic(U);
end

function E=struct_raw(W,phi,psi,internal)
    % Raw structured word with the complete shifted pairing: D'_k = G_k * Q_k, Q_k = product of internal
    % shifted sections on planes (2,3),(4,5),...,(M-2,M-1); internal: (M/2-1) x K angles.
    M=size(W{1},1);K=numel(W)-1;U={W{1}};
    for k=1:K
        Q=eye(M);for r=1:M/2-1,Q=Q*rot(M,2*r,2*r+1,internal(r,k));end
        U{end+1}=rot(M,1,2,phi(k));U{end+1}=rot(M,1,M,-psi(k))*Q;U{end+1}=W{k+1}; %#ok<AGROW>
    end
    E=basic(U);
end

function x=partial_x(U)
    % x_j = P_j' e_1 with P_j = U_j ... U_0, j = 0..L-1 (basic reduced sphere coordinates).
    M=size(U{1},1);P=eye(M);x=zeros(M,numel(U)-1);
    for j=1:numel(U)-1,P=U{j}*P;x(:,j)=P'*[1;zeros(M-1,1)];end
end

function F=gram(x)
    % Layerwise Gram function F_k = <a,c> - <a,b><b,c> with (a,b,c) = (x_{3k-3},x_{3k-2},x_{3k-1}).
    K=size(x,2)/3;F=zeros(1,K);
    for k=1:K,a=x(:,3*k-2);b=x(:,3*k-1);c=x(:,3*k);F(k)=a'*c-(a'*b)*(b'*c);end
end

function h=filters(E)
    % Type-I analysis filters: h_m[M*l+p] = E_mp[l]; returns M x M(L+1).
    [M,~,L1]=size(E);h=zeros(M,M*L1);
    for l=0:L1-1,for p=0:M-1,h(:,M*l+p+1)=E(:,p+1,l+1);end,end
end

function r=paraunitarity(E)
    % || sum_l E_l' E_{l+s} - delta_s I || over all lags s (exact polynomial paraunitarity check).
    [M,~,L1]=size(E);r=0;
    for s=-(L1-1):(L1-1)
        S=zeros(M);for l=1:L1,j=l+s;if j>=1&&j<=L1,S=S+E(:,:,l)'*E(:,:,j);end,end
        r=max(r,norm(S-(s==0)*eye(M),'fro'));
    end
end

function A=block_matrix(E,N)
    % Periodic finite-block analysis operator of length N (N divisible by M), correlation convention
    % y_m[k] = sum_n h_m[n] x[(M*k+n) mod N], outputs ordered subband-major (index m*B+k).
    % Rows of A are the analysis filters h_m (forward) at offsets M*k; U = A' has the synthesis basis
    % functions h_m as columns, matching the current fir_mallat_r22a convention (columns = conv(pulse,h)).
    M=size(E,1);assert(mod(N,M)==0);B=N/M;h=filters(E);A=zeros(N);
    for m=1:M,for k=0:B-1,for n=0:size(h,2)-1,c=mod(M*k+n,N)+1;A((m-1)*B+k+1,c)=A((m-1)*B+k+1,c)+h(m,n+1);end,end,end
end

function y=apply_polyphase(E,x)
    % Block convolution: y_blk[k] = sum_l E_l x_blk[(k+l) mod B], x_blk[k]_p = x[M*k+p].
    [M,~,L1]=size(E);N=numel(x);B=N/M;X=reshape(x,M,B);Y=zeros(M,B);
    for k=0:B-1,for l=0:L1-1,Y(:,k+1)=Y(:,k+1)+E(:,:,l+1)*X(:,mod(k+l,B)+1);end,end
    y=reshape(Y.',[],1);
end

function y=apply_filterbank(E,x)
    % Circular cross-correlation with each analysis filter, then decimation by M.
    M=size(E,1);N=numel(x);B=N/M;h=filters(E);y=zeros(N,1);
    for m=1:M
        full=zeros(N,1);for t=0:N-1,for n=0:size(h,2)-1,full(t+1)=full(t+1)+h(m,n+1)*x(mod(t+n,N)+1);end,end
        y((m-1)*B+(1:B))=full(1:M:end);
    end
end

function Wt=w1_terminal(R3,T)
    % Terminal factor enforcing W1 for E(1) = Wt*T: Wt*T*v = e1. Wt = diag(1,R3)*A with A in SO(M), A*(T*v)=e1.
    M=size(T,1);u=T*ones(M,1)/sqrt(M);e=[1;zeros(M-1,1)];
    if norm(u-e)<1e-14,A=eye(M);else
        w=(u-e)/norm(u-e);A=eye(M)-2*(w*w');            % reflection u->e1 (det -1)
        A=diag([1,-ones(1,1),ones(1,M-2)])*A;           % compose with a reflection fixing e1 -> det +1, still maps u->e1
    end
    Wt=blkdiag(1,R3)*A;
end
