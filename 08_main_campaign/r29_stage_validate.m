function result=r29_stage_validate(paths,p)
% VALIDATION selection and candidate freeze (sec. 17). Uses DESIGN (support learning) and VALIDATION only.
% Per family: the 3 best distinct DESIGN finals are scored on VALIDATION with their DESIGN supports;
% the minimum VALIDATION J_broad is frozen as the primary candidate. s_ref (sec. 14) is a diagnostic only.
% HYBRID-M4-DYADIC reuses the frozen M4-INDEPENDENT-STRUCT-(1,1) parameters (no separate optimization).
    out=fullfile(paths.c,'data');freeze_file=fullfile(paths.c,'config','r29_candidate_freeze.mat');
    assert(~isfile(freeze_file),'candidate freeze already exists; never overwritten');
    Gd=r29_io('load_split',paths,'design');Gv=r29_io('load_split',paths,'validation');
    controls=r29_io('controls',paths);ref=controls(strcmp({controls.name},'INDEPENDENT_PUFB'));
    eref=split_error(ref.U,[],Gd,Gv,p.rho_design);rows={};cand=struct('name',{},'label',{},'params',{},'Unu',{},'Utau',{},'rho',{});
    for i=1:numel(p.families_all)
        name=p.families_all{i};f=r29_families('spec',name);
        if strcmp(name,'HYBRID-M4-DYADIC')
            base=cand(strcmp({cand.name},'M4-INDEPENDENT-STRUCT-(1,1)')&strcmp({cand.label},'BROAD'));P=base.params;Jd=nan;src="FROZEN_(1,1)";
        else
            [P,Jd]=distinct(union_runs(out,['r29c_broad_',r29_io('tag',name),'.mat'],['r29c_refine_',r29_io('tag',name),'.mat']),3);src="DESIGN_TOP3_UNION";
        end
        scores=zeros(size(P,1),1);sref=scores;
        for j=1:size(P,1)
            [~,ax]=r29_families('U',f,P(j,:));e=split_error(ax{1},ax{2},Gd,Gv,p.rho_design);
            scores(j)=mean(log10(max(e,1e-16)));sref(j)=mean(e./(e+eref));
            rows(end+1,:)={string(name),"BROAD",src,j,Jd(min(j,end)),scores(j),sref(j),e}; %#ok<AGROW>
        end
        [~,b]=min(scores);[~,bs]=min(sref);[~,ax]=r29_families('U',f,P(b,:));
        cand(end+1)=struct('name',name,'label','BROAD','params',P(b,:),'Unu',ax{1},'Utau',ax{2},'rho',p.rho_design); %#ok<AGROW>
        sens(i)=struct('family',name,'selected_by_J',b,'selected_by_sref',bs,'same',b==bs); %#ok<AGROW>
    end
    % Diagnostic budget-specific candidates: best VALIDATION e(rho) among their DESIGN finals.
    for name=["M4-INDEPENDENT-STRUCT-(1,1)","M4-BASIC-(3,3)"]
        f=r29_families('spec',char(name));
        for rho=p.rho_diag
            tg=sprintf('%s_rho%03d.mat',r29_io('tag',char(name)),round(100*rho));[P,Jd]=distinct(union_runs(out,['r29c_budget_',tg],['r29c_refine_budget_',tg]),3);
            ev=zeros(size(P,1),1);for j=1:size(P,1),[~,ax]=r29_families('U',f,P(j,:));ev(j)=split_error(ax{1},ax{2},Gd,Gv,rho);
                rows(end+1,:)={name,sprintf("DIAG_RHO_%.2f",rho),"DESIGN_TOP3",j,Jd(j),log10(ev(j)),nan,nan(1,5)};end %#ok<AGROW>
            [~,b]=min(ev);[~,ax]=r29_families('U',f,P(b,:));
            cand(end+1)=struct('name',char(name),'label',sprintf('DIAG_RHO_%.2f',rho),'params',P(b,:),'Unu',ax{1},'Utau',ax{2},'rho',rho); %#ok<AGROW>
        end
    end
    % Historical/control context on VALIDATION (not selected; DESIGN supports).
    for c=controls
        e=split_error(c.U,[],Gd,Gv,p.rho_design);rows(end+1,:)={string(c.name),"CONTROL","FIXED",1,nan,mean(log10(max(e,1e-16))),mean(e./(e+eref)),e}; %#ok<AGROW>
    end
    T=cell2table(rows,'VariableNames',{'Family','Label','Source','Candidate','JDesign','JValidation','SRefValidation','eValidation'});
    writetable(T,fullfile(out,'r29c_validation_scores.csv'));writetable(struct2table(sens),fullfile(out,'r29c_objective_sensitivity.csv'));
    gates=predeclared_gates(p);
    save(freeze_file,'cand','gates');
    man=struct('frozen_utc',char(datetime('now','TimeZone','UTC')),'freeze_mat_sha256',r22a_io('hash',freeze_file), ...
        'test_opened',false,'support_learning','mean |U^H G U|^2 over the 300 DESIGN channels, stable descending sort, top round(rho*n^2)', ...
        'tree','per family spec (r29_families.m); HYBRID uses fir_mallat_r22a(6,pi/4) on the 6-sample M_tau lowpass block', ...
        'candidates',struct('name',{cand.name},'label',{cand.label},'params',{cand.params},'Unu_sha',cellfun(@(x)r22dcb_math('arrayhash',x),{cand.Unu},'UniformOutput',false), ...
        'Utau_sha',cellfun(@(x)r22dcb_math('arrayhash',x),{cand.Utau},'UniformOutput',false)),'gates',gates, ...
        'sources_sha256',struct('r29_families',r22a_io('hash',which('r29_families')),'r29_pu',r22a_io('hash',which('r29_pu')),'r29_stage_design',r22a_io('hash',which('r29_stage_design')),'r29_stage_refine',r22a_io('hash',which('r29_stage_refine'))));
    r22a_io('json',fullfile(paths.c,'config','r29_candidate_freeze_manifest.json'),man);
    result=struct('candidates',numel(cand),'objective_sensitivity_same',all([sens.same]));disp(T(:,[1 2 4 5 6 7]));
end

function runs=union_runs(out,varargin)
    % Union of original (fminunc) and protocol-revision (pattern search) DESIGN finals.
    J=[];P=[];for i=1:numel(varargin),l=load(fullfile(out,varargin{i}),'runs');J=[J;l.runs.JFinal];P=[P;l.runs.Params];end %#ok<AGROW>
    runs=table(J,P,'VariableNames',{'JFinal','Params'});
end

function [P,J]=distinct(runs,k)
    [J,o]=sort(runs.JFinal);P=runs.Params(o,:);keep=true(size(J));
    for i=2:numel(J),for j=1:i-1,if keep(j)&&(abs(J(i)-J(j))<1e-7||norm(mod(P(i,:)-P(j,:)+pi,2*pi)-pi)<1e-4),keep(i)=false;end,end,end
    P=P(keep,:);J=J(keep);P=P(1:min(k,end),:);J=J(1:min(k,end));
end

function e=split_error(A,B,Gs,Ge,rhos)
    % Support from the support bank Gs (DESIGN), normalized omitted energy averaged over the evaluation bank Ge.
    sup=r29_io('mean_energy',A,B,Gs);[~,order]=sort(sup(:),'descend');
    nrm=r29_io('mean_normalized_energy',A,B,Ge);cum=cumsum(nrm(order));e=max(1-cum(round(rhos*288^2)).',0);
end

function g=predeclared_gates(p)
    % Operational thresholds fixed BEFORE TEST is opened (sec. 23-24 left them qualitative).
    g=struct();
    g.comparator_set={'IDENTITY','INDEPENDENT_PUFB','SHARED_PUFB','HAAR_X_HAAR','D4_X_D4','RECT_OTFS','DYADIC-W1-D1','DYADIC-W1-D3','DYADIC-W1-D7'};
    g.comparator_rule='at each rho_RX, the comparator with the lowest TEST mean of the metric (re-chosen per rho and per metric)';
    g.tested_architectures={'M4-INDEPENDENT-STRUCT-(1,1)','M4-BASIC-(3,3)','M4-BASIC-D1','M4-INDEPENDENT-STRUCT-(1,2)','M4-BASIC-(3,6)','M4-SHARED-STRUCT-(1,1)'};
    g.diagnostic_only={'HYBRID-M4-DYADIC','DIAG_RHO_* budget-specific candidates'};
    g.practical_grid=[0.30,0.35,0.40,0.45,0.50,0.60,0.70];
    g.uncertainty='paired channel bootstrap within velocity, 10000 replicates, seed 292200000, percentile 95% CI of the relative reduction (ratio of means)';
    g.material='relative reduction >= 10% with the 95% CI excluding 0';
    g.ber_gate_A='DeltaM (arm R, regime-balanced 3-24 dB): material at two adjacent practical_grid points for one tested architecture';
    g.ber_gate_B=['support regret R_E: relative reduction >= 25% (CI excludes 0) AND eta_capture (pooled) increase with CI excluding 0, ', ...
        'both at the same two adjacent practical_grid points, AND the R_E reduction has the same sign in >= 4 of 5 velocities at both points'];
    g.classes=struct('M4_A','no material improvement at any practical_grid point in oracle omitted energy, offline omitted energy, R_E, eta_capture, Jaccard, or DeltaM', ...
        'M4_B','oracle omitted energy material at >= 2 adjacent points but offline omitted energy, eta_capture and R_E not material at any 2 adjacent points', ...
        'M4_C','offline omitted energy, eta_capture or R_E material at >= 3 adjacent practical_grid points', ...
        'M4_D','material only at 1-2 adjacent points or only outside the practical grid, or the best-representation crossover rho moves by >= one grid step');
    g.rho_test=p.rho_test;
end
