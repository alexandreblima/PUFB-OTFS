function varargout=r29_io(action,varargin)
% Shared I/O and energy-map helpers for R2.9 VALIDATE/TEST stages.
    switch action
        case 'load_split',varargout{1}=load_split(varargin{:});
        case 'tag',varargout{1}=regexprep(varargin{1},'[^A-Za-z0-9]+','_');
        case 'controls',varargout{1}=controls(varargin{:});
        case 'coeffs',varargout{1}=coeffs(varargin{:});
        case 'mean_energy',a=energy_maps(varargin{:});varargout{1}=a{1};
        case 'mean_normalized_energy',a=energy_maps(varargin{:});varargout{1}=a{2};
        otherwise,error('R29IO:Action','Unknown action %s.',action);
    end
end

function G=load_split(paths,split)
    file=fullfile(paths.data,sprintf('r29_%s_manifest.csv',split));
    T=readtable(file,'TextType','string','Delimiter',',','VariableNamingRule','preserve');
    seal=jsondecode(fileread(fullfile(paths.data,sprintf('r29_%s_seal.json',split))));
    assert(strcmp(seal.manifest_sha256,r22a_io('hash',file)),'manifest changed');
    if split=="test",assert(isfile(fullfile(paths.c,'config','r29_candidate_freeze_manifest.json')),'TEST cannot be opened before the candidate freeze');end
    G=complex(zeros(288,288,height(T)));
    for i=1:height(T),l=load(T.AbsolutePath(i),'G');G(:,:,i)=l.G;end
end

function C=controls(paths)
    % Fixed historical controls, identical matrices to R2.6 (asserted there against R2.2D/R2.2E/R2.2H).
    r22e=fullfile(paths.v2,'results','r2_2e_receiver_floor');r22h=fullfile(paths.v2,'results','r2_2h_natural_baselines');
    frozen=load(fullfile(r22e,'r22e_supports.mat'),'methods');frozen=frozen.methods;h=load(fullfile(r22h,'r22h_frozen_transforms.mat'),'U');
    assert(norm(frozen(3).U-h.U{5},'fro')<1e-13);
    names={'IDENTITY','INDEPENDENT_PUFB','SHARED_PUFB','HAAR_X_HAAR','D4_X_D4','RECT_OTFS'};U={eye(288),frozen(3).U,frozen(2).U,h.U{4},h.U{3},h.U{2}};
    C=struct('name',names,'U',U);
end

function H=coeffs(A,B,G)
    % A,B = (Unu,Utau) for separable real families; B = [] means A is a dense 288x288 U.
    if isempty(B)
        H=pagemtimes(pagemtimes(A,'ctranspose',G,'none'),A);
    else
        H=r29_families('transform',G,A,B);
    end
end

function maps=energy_maps(A,B,G)
    C=size(G,3);sup=zeros(288);nrm=zeros(288);
    for s=1:50:C
        idx=s:min(C,s+49);a=abs(coeffs(A,B,G(:,:,idx))).^2;sup=sup+sum(a,3);nrm=nrm+sum(a./sum(a,[1,2]),3);
    end
    maps={sup/C,nrm/C};
end
