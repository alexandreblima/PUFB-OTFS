function result=run_r29(root,stage,varargin)
% R2.9 continuation driver. Stages: CONTROLS, DIMAUDIT, DATA, DESIGN, VALIDATE, TEST, PAPR.
% TEST channels are never used before the candidate-freeze manifest exists.
    [paths,p]=setup(root);
    switch stage
        case 'CONTROLS',result=controls(paths,p);
        case 'DIMAUDIT',result=dimaudit(paths,p);
        case 'DATA',result=make_data(paths,p);
        otherwise
            fn=str2func(['r29_stage_',lower(stage)]);result=fn(paths,p,varargin{:});
    end
end

function [paths,p]=setup(root)
    v2=fullfile(root,'pufb_otfs_v2');base=fullfile(v2,'results','next_program','phase1f_multichannel_pufb');
    pkg=fileparts(mfilename('fullpath'));pkg=fileparts(pkg);
    paths=struct('root',root,'v2',v2,'audit',fullfile(pkg,'08_main_campaign'),'base',base, ...
        'b',fullfile(base,'r2_9b_m4_implementation'),'c',fullfile(base,'r2_9c_design'),'d',fullfile(base,'r2_9d_compression_support'), ...
        'f',fullfile(base,'r2_9f_papr'),'h',fullfile(base,'r2_9h_scientific_closure'),'data',fullfile(base,'channels'));
    addpath(paths.audit,fullfile(pkg,'02_shared_lib'),fullfile(pkg,'03_axis_specific_pufb'),fullfile(pkg,'04_ber_confirmation'), ...
        fullfile(pkg,'06_receiver_floor'),fullfile(pkg,'01_channel_generator'),fullfile(v2,'audit','r2_5_abid_full_mmse'));
    for d={paths.c,paths.d,paths.f,paths.h},for s={'','data','figures','scripts','config'},if ~isfolder(fullfile(d{1},s{1})),mkdir(fullfile(d{1},s{1}));end,end,end
    if ~isfolder(paths.data),mkdir(paths.data);end
    p=jsondecode(fileread(fullfile(base,'r29_seed_plan.json')));
    p.rho_design=[0.30,0.40,0.50,0.60,0.70];p.rho_test=[0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.00];
    p.rho_diag=[0.40,0.50,0.70];
    p.families={'M4-INDEPENDENT-STRUCT-(1,1)','M4-BASIC-(3,3)','DYADIC-W1-D3','DYADIC-W1-D7','M4-BASIC-D1', ...
        'M4-INDEPENDENT-STRUCT-(1,2)','M4-BASIC-(3,6)','M4-SHARED-STRUCT-(1,1)'};   % optimization priority order (sec. 29)
    p.families_all=[p.families,{'DYADIC-W1-D1','HYBRID-M4-DYADIC'}];   % D1 = re-optimized current-architecture control; HYBRID not optimized
end

%% ------------------------------------------------------------------ CONTROLS
function result=controls(paths,~)
    rng(293001,'twister');rows=[];
    for L=[1,3,7]
        f=r29_families('spec',sprintf('DYADIC-W1-D%d',L));
        for t=1:20
            c=2*pi*rand(1,L);E=r29_families('axis_E',f.axes{1},c);detp=det_poly(E);
            A12=r29_families('axis_U',f.axes{1},c,12);A24=r29_families('axis_U',f.axes{2},c,24);x=randn(24,1);
            single=r29_pu('block_matrix',E,24);
            rows=[rows;L,r29_pu('paraunitarity',E),norm(r29_pu('E1',E)*ones(2,1)/sqrt(2)-[1;0]),norm(A12'*A12-eye(12),'fro'),norm(A24'*A24-eye(24),'fro'), ...
                max(abs(single*x-r29_pu('apply_polyphase',E,x))),max(abs(single*x-r29_pu('apply_filterbank',E,x))),abs(norm(A24'*x)^2-norm(x)^2)/norm(x)^2, ...
                size(E,3)-1,nnz(abs(r29_pu('filters',E))>1e-13)/2,mcmillan(detp,L)]; %#ok<AGROW>
        end
    end
    T=array2table(rows,'VariableNames',{'Depth','Paraunitarity','W1','TreeOrth12','TreeOrth24','PolyphaseVsMatrix','FilterbankVsMatrix','EnergyRel', ...
        'PolyOrder','NonzeroTapsPerFilter','McMillanResidual'});
    % D1 must reproduce the current PUFB transform up to coefficient signs (orientation convention).
    sign_err=0;f=r29_families('spec','DYADIC-W1-D1');
    for alpha=[0,pi/4,pi/3,0.78797323934301,0.4,1.3,2.2,2.9]
        [~,~,meta]=fir_mallat_r22a(6,alpha);E=zeros(2,2,2);for l=0:1,E(:,:,l+1)=[meta.h(2*l+1),meta.h(2*l+2);meta.g(2*l+1),meta.g(2*l+2)];end
        x0=null(E(:,:,2));theta=atan2(x0(2),x0(1));
        [P12,~,~]=fir_mallat_r22a(12,alpha);[P24,~,~]=fir_mallat_r22a(24,alpha);
        sign_err=max([sign_err,colsign(r29_families('axis_U',f.axes{1},theta,12),P12),colsign(r29_families('axis_U',f.axes{2},theta,24),P24)]);
    end
    % Fast separable transform equals the dense product.
    G=randn(288,288,3)+1i*randn(288,288,3);f=r29_families('spec','M4-INDEPENDENT-STRUCT-(1,2)');c=randn(1,r29_families('nparams',f));
    [U,ax]=r29_families('U',f,c);H=r29_families('transform',G,ax{1},ax{2});tr=0;for i=1:3,tr=max(tr,max(abs(H(:,:,i)-U'*G(:,:,i)*U),[],'all'));end
    summary=struct('max',varfun(@max,T(:,2:8)),'support_by_depth',unique(T(:,{'Depth','PolyOrder','NonzeroTapsPerFilter'})), ...
        'D1_equals_current_up_to_signs',sign_err,'fast_transform_vs_dense',tr);
    writetable(T,fullfile(paths.b,'data','r29b_dyadic_controls.csv'));
    result=struct('max_paraunitarity',max(T.Paraunitarity),'max_w1',max(T.W1),'max_tree_orth',max([T.TreeOrth12;T.TreeOrth24]), ...
        'max_polyphase_matrix',max(T.PolyphaseVsMatrix),'max_filterbank_matrix',max(T.FilterbankVsMatrix),'max_energy',max(T.EnergyRel), ...
        'max_mcmillan_residual',max(T.McMillanResidual),'support_depth3',2*(3+1),'support_depth7',2*(7+1), ...
        'taps_depth3',max(T.NonzeroTapsPerFilter(T.Depth==3)),'taps_depth7',max(T.NonzeroTapsPerFilter(T.Depth==7)), ...
        'D1_equals_current_up_to_signs',sign_err,'fast_transform_vs_dense',tr);
    result.PASS=result.max_paraunitarity<1e-12&&result.max_w1<1e-12&&result.max_tree_orth<1e-12&&result.max_polyphase_matrix<1e-12 ...
        &&result.max_filterbank_matrix<1e-12&&result.max_energy<1e-12&&result.max_mcmillan_residual<1e-10&&sign_err<1e-12&&tr<1e-10;
    r22a_io('json',fullfile(paths.b,'data','r29b_dyadic_controls_gate.json'),result);clear summary
end

function r=colsign(A,B)
    % Columnwise equality up to a sign flip per basis function.
    r=max(min(vecnorm(A-B),vecnorm(A+B)));
end

function r=mcmillan(d,L)
    % det E = +-z^-L exactly (McMillan degree L): one unit coefficient at lag L, zeros elsewhere.
    r=abs(abs(d(L+1))-1)+max(abs(d([1:L,L+2:end])),[],'all');if isempty(r),r=0;end
end

function d=det_poly(E)
    % det of a 2x2 polynomial matrix as coefficient vector in z^-1 (McMillan degree check: det = +-z^-L).
    a=squeeze(E(1,1,:)).';b=squeeze(E(1,2,:)).';c=squeeze(E(2,1,:)).';e=squeeze(E(2,2,:)).';d=conv(a,e)-conv(b,c);
end

%% ------------------------------------------------------------------ DIMAUDIT
function result=dimaudit(paths,~)
    rng(293002,'twister');cases={ ...
        'M4S',1,'SINGLE',12;'M4S',1,'SINGLE',24;'M4S',2,'SINGLE',12;'M4S',2,'SINGLE',24; ...
        'M4B',1,'SINGLE',12;'M4B',1,'SINGLE',24;'M4B',3,'SINGLE',12;'M4B',3,'SINGLE',24;'M4B',6,'SINGLE',12;'M4B',6,'SINGLE',24; ...
        'DY',1,'MALLAT',12;'DY',1,'MALLAT',24;'DY',3,'MALLAT',12;'DY',3,'MALLAT',24;'DY',7,'MALLAT',12;'DY',7,'MALLAT',24;'M4S',1,'HYBRID',24};
    thresholds=[1e-6,1e-8,1e-10,1e-12];rows=[];
    for i=1:size(cases,1)
        a=struct('kind',cases{i,1},'depth',cases{i,2},'tree',cases{i,3});N=cases{i,4};
        switch a.kind,case 'DY',d=a.depth;case 'M4B',d=3*a.depth+3;case 'M4S',d=8*a.depth+3;end
        for point=1:5
            c=randn(1,d);if a.kind=="M4S",K=a.depth;c(6*K+(1:2*K))=0.3+(pi-0.6)*rand(1,2*K);end
            if a.kind=="M4B",for j=1:a.depth,c(3*j-2:3*j-1)=0.3+(pi-0.6)*rand(1,2);end,end   % away from hyperspherical poles
            fu=@(x)reshape(r29_families('axis_U',a,x,N),[],1);[s,J]=jac_sv(fu,c,1e-3);[~,J2]=jac_sv(fu,c,5e-4);fd=norm(J-J2)/norm(J);
            rows=[rows;{a.kind,a.depth,a.tree,N,point,d},num2cell(arrayfun(@(t)nnz(s>t*s(1)),thresholds)),{min(s)/s(1),fd}]; %#ok<AGROW>
        end
    end
    T=cell2table(rows,'VariableNames',{'Kind','Depth','Tree','N','Point','DAbstract','DFinite_1e6','DFinite_1e8','DFinite_1e10','DFinite_1e12','MinRelSV','FDRelError'});
    T.RFinite_1e8=T.DFinite_1e8./T.DAbstract;writetable(T,fullfile(paths.b,'data','r29b_finite_dimension_audit.csv'));
    result=groupsummary(T,{'Kind','Depth','Tree','N','DAbstract'},{'min','max'},{'DFinite_1e6','DFinite_1e8','DFinite_1e10','DFinite_1e12','MinRelSV','FDRelError'});disp(result);
    writetable(result,fullfile(paths.b,'data','r29b_finite_dimension_summary.csv'));
end

function [s,J]=jac_sv(f,c,h)
    % Five-point central differences (truncation ~h^4) for an accurate Jacobian spectrum.
    f0=f(c);J=zeros(numel(f0),numel(c));
    for i=1:numel(c)
        e=zeros(size(c));e(i)=h;J(:,i)=(-f(c+2*e)+8*f(c+e)-8*f(c-e)+f(c-2*e))/(12*h);
    end
    s=svd(J);
end

%% ------------------------------------------------------------------ DATA
function result=make_data(paths,p)
    physical=jsondecode(fileread(fullfile(paths.v2,'audit','r2_0a_abid_true2d_pilot','r20a_preregistration.json')));
    splits={'DESIGN','VALIDATION','TEST'};rows=cell(0,6);
    for s=1:3
        split=splits{s};directory=fullfile(paths.data,lower(split));if ~isfolder(directory),mkdir(directory);end
        base=p.bases.(split);count=p.per_velocity.(split);
        for v=1:5
            for c=1:count
                seed=base+1000*v+c;path=fullfile(directory,sprintf('v%d_c%03d.mat',v,c));
                if ~isfile(path),G=generate_channel_r20a(24,12,p.velocities_mps(v),seed,physical);save(path,'G','seed','-v7');end
                rows(end+1,:)={string(split),v,c,seed,string(path),string(r22a_io('hash',path))}; %#ok<AGROW>
            end
        end
        fprintf('R29 DATA %s generated\n',split);
    end
    T=cell2table(rows,'VariableNames',{'Split','VelocityIndex','ChannelID','Seed','AbsolutePath','SHA256'});
    files=dir(fullfile(paths.v2,'results','**','*manifest*.csv'));overlap=0;
    for i=1:numel(files)
        if contains(files(i).folder,'phase1f_multichannel_pufb'),continue;end
        M=readtable(fullfile(files(i).folder,files(i).name),'Delimiter',',','VariableNamingRule','preserve');
        if ismember('Seed',M.Properties.VariableNames),overlap=overlap+numel(intersect(T.Seed,M.Seed));end
    end
    for s=1:3
        S=T(T.Split==splits{s},:);path=fullfile(paths.data,sprintf('r29_%s_manifest.csv',lower(splits{s})));writetable(S,path);
        r22a_io('json',fullfile(paths.data,sprintf('r29_%s_seal.json',lower(splits{s}))),struct('split',splits{s},'channels',height(S), ...
            'manifest_sha256',r22a_io('hash',path),'sealed_utc',char(datetime('now','TimeZone','UTC')),'opened_for_metrics',false));
    end
    result=struct('channels',height(T),'seed_overlap_with_historical',overlap,'generator','generate_channel_r20a (unchanged)');
    r22a_io('json',fullfile(paths.data,'r29_data_gate.json'),result);assert(overlap==0);
end
