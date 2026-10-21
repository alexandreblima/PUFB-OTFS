function result=r29_stage_ber(paths,p,mode)
% Gated BER campaign (sec. 25). Runs only if r29d_gates_classification.json reports BER_GATE_PASS.
% Paired channels (R2.9 TEST, 1000), bits, symbols, noise and frame counts across all representations and budgets.
% Receiver/detector: R2.6/R2.5 kernel copied verbatim (mismatched LMMSE, lambda = sigma^2, mrg32k3a substreams per 16-frame batch).
    out=fullfile(paths.base,'r2_9e_optional_ber');for s={'','data','figures','scripts','config'},if ~isfolder(fullfile(out,s{1})),mkdir(fullfile(out,s{1}));end,end
    gate=jsondecode(fileread(fullfile(paths.d,'data','r29d_gates_classification.json')));assert(gate.BER_GATE_PASS,'BER gate did not pass');
    an=jsondecode(fileread(fullfile(paths.d,'data','r29d_analytic_gate.json')));assert(an.PASS);
    switch mode
        case 'DECLARE'
            file=fullfile(out,'config','r29e_ber_protocol.json');assert(~isfile(file));
            r22a_io('json',file,struct('declared_utc',char(datetime('now','TimeZone','UTC')), ...
                'rho',[0.40,0.50,0.60,0.70,1.00],'midpoint',0.60, ...
                'midpoint_justification',['Gate A was met at rho = 0.50 and 0.60 (adjacent practical-grid points) by M4-BASIC-(3,3) and ', ...
                    'M4-INDEPENDENT-STRUCT-(1,2); 0.60 is where the analytic DeltaM advantage is resolved but is not in the base set {0.4,0.5,0.7,1.0}. Declared before any BER.'], ...
                'representations',{{'IDENTITY','INDEPENDENT_PUFB','RECT_OTFS','DYADIC-W1-D1','DYADIC-W1-D7','M4-BASIC-(3,3)','M4-INDEPENDENT-STRUCT-(1,2)','M4-INDEPENDENT-STRUCT-(1,1)'}}, ...
                'snr_dB',0:2:30,'frames_per_channel_per_snr',32,'batches_of_16',2, ...
                'seeds','data 294100000 + 100000*v + 1000*snr + channel; noise 294200000 + 100000*v + 1000*snr + channel (paired across methods)', ...
                'supports','frozen DESIGN offline orders (r29d_offline_orders.mat)','target_BER',1e-3, ...
                'statistics','r25_stats.summarize (first descending bracket crossing, no extrapolation, paired channel bootstrap within velocity)'));
            result=true;
        case 'RUN',result=run(paths,p,out);
    end
end

function result=run(paths,p,out)
    proto=jsondecode(fileread(fullfile(out,'config','r29e_ber_protocol.json')));
    o=load(fullfile(paths.d,'data','r29d_offline_orders.mat'));l=load(fullfile(paths.c,'config','r29_candidate_freeze.mat'),'cand');C=r29_io('controls',paths);
    methods=struct('Name',{},'U',{},'Mask',{});
    for k=1:numel(proto.rho)
        K=round(proto.rho(k)*288^2);
        for name=string(proto.representations)'
            j=find(string(o.names)==name);assert(isscalar(j));
            c=find(strcmp({C.name},name));if ~isempty(c),U=C(c).U;else,f=find(strcmp({l.cand.name},name)&strcmp({l.cand.label},'BROAD'));U=kron(l.cand(f).Unu,l.cand(f).Utau);end
            mask=false(288);mask(o.order(1:K,j))=true;methods(end+1)=struct('Name',sprintf('%s_RHO%03d',name,round(100*proto.rho(k))),'U',U,'Mask',mask); %#ok<AGROW>
        end
    end
    T=readtable(fullfile(paths.data,'r29_test_manifest.csv'),'TextType','string','Delimiter',',','VariableNamingRule','preserve');
    directory=fullfile(out,'data','ber_points');if ~isfolder(directory),mkdir(directory);end
    if isempty(gcp('nocreate')),parpool('Processes',6);end
    for v=1:5
        files=T.AbsolutePath(T.VelocityIndex==v);n=numel(files);
        for snr=proto.snr_dB(:)'
            path=fullfile(directory,sprintf('v%d_snr%02d.mat',v,snr));if isfile(path),continue;end
            variance=10^(-snr/10);local=(1:n)';ds=294100000+100000*v+1000*snr+local;ns=294200000+100000*v+1000*snr+local;errors=zeros(numel(methods),n);
            parfor ch=1:n
                g=load(files(ch),'G');errors(:,ch)=sparse_errors(g.G,methods,variance,ds(ch),ns(ch),proto.batches_of_16); %#ok<PFBNS>
            end
            bits=576*16*proto.batches_of_16;names={methods.Name};save(path,'errors','bits','names');fprintf('R29 BER v%d SNR%02d\n',v,snr);
        end
    end
    result=numel(methods);
end

function errors=sparse_errors(G,methods,variance,data_seed,noise_seed,batches)
    % Verbatim R2.6 kernel (same bit/symbol/noise generation as R2.5: mrg32k3a, substream per 16-frame batch).
    count=numel(methods);signal=cell(count,1);noise=cell(count,1);
    for index=1:count
        H=methods(index).U'*G*methods(index).U;Hs=H.*methods(index).Mask;
        factor=chol(Hs'*Hs+variance*eye(288));W=factor\(factor'\Hs');signal{index}=W*H;noise{index}=W*methods(index).U';
    end
    data=RandStream('mrg32k3a','Seed',data_seed);stream=RandStream('mrg32k3a','Seed',noise_seed);errors=zeros(count,1);
    for batch=1:batches
        data.Substream=batch;stream.Substream=batch;
        bits=logical(randi(data,[0,1],576,16));real_bits=bits(1:2:end,:);imag_bits=bits(2:2:end,:);
        symbols=(1-2*double(real_bits)+1i*(1-2*double(imag_bits)))/sqrt(2);
        physical=sqrt(variance/2)*(randn(stream,288,16)+1i*randn(stream,288,16));
        for index=1:count
            estimate=signal{index}*symbols+noise{index}*physical;
            errors(index)=errors(index)+nnz((real(estimate)<0)~=real_bits)+nnz((imag(estimate)<0)~=imag_bits);
        end
    end
end
