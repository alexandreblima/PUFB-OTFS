function result=r22dhma_kernel(H,mask,snr_db,p,reference_check)
    dimension=size(H,1);identity=eye(dimension);Hs=H.*mask;omitted=H-Hs;
    assert(norm(H-Hs-omitted,'fro')<=p.identity_tolerance*(1+norm(H,'fro')));
    [left,diagonal,right]=svd(Hs,'econ');singular=diag(diagonal);omitted_singular=svd(omitted);
    maximum=singular(1);minimum=singular(end);tolerance=p.rank_relative_tolerance*maximum;active=singular>tolerance;
    condition=Inf;if minimum>0,condition=maximum/minimum;end
    pseudoinverse=right(:,active)*(left(:,active)'./singular(active));
    limit_error=pseudoinverse*H-identity;limit_omitted=pseudoinverse*omitted;
    result=struct('Singular',singular,'Static',[maximum,minimum,condition,nnz(active),tolerance, ...
        norm(omitted,'fro')/norm(H,'fro'),norm(omitted,'fro'),omitted_singular(1)], ...
        'Limit',[norm(limit_error,'fro')^2/dimension,norm(limit_omitted,'fro')^2/dimension], ...
        'Values',zeros(numel(snr_db),23));
    sparse_gram=Hs'*Hs;full_gram=H'*H;
    for point=1:numel(snr_db)
        variance=10^(-snr_db(point)/10);system=sparse_gram+variance*identity;factor=chol(system);
        receiver=factor\(factor'\Hs');regularization=receiver*Hs-identity;amplified=receiver*omitted;
        reg=norm(regularization,'fro')^2/dimension;omit=norm(amplified,'fro')^2/dimension;
        cross=2*real(sum(regularization(:).*conj(amplified(:))))/dimension;
        receiver_f=norm(receiver,'fro');noise=variance*receiver_f^2/dimension;
        direct=(norm(receiver*H-identity,'fro')^2+variance*receiver_f^2)/dimension;
        consistent=reg+noise;penalty=direct-consistent;reconstruction=abs(direct-reg-omit-cross-noise);
        penalty_error=abs(penalty-omit-cross);
        assert(max(reconstruction,penalty_error)<=p.identity_tolerance*(1+abs(direct)),'R22DHMA:Decomposition','Exact decomposition failed.');
        full_factor=chol(full_gram+variance*identity);full_receiver=full_factor\(full_factor'\H');
        full_mse=(norm(full_receiver*H-identity,'fro')^2+variance*norm(full_receiver,'fro')^2)/dimension;
        amplified_f=norm(amplified,'fro');amplified_2=norm(amplified,2);gain_f=NaN;gain_2=NaN;
        if norm(omitted,'fro')>0,gain_f=amplified_f/norm(omitted,'fro');end
        if omitted_singular(1)>0,gain_2=amplified_2/omitted_singular(1);end
        lambda_min=minimum^2+variance;lambda_max=maximum^2+variance;condition_A=lambda_max/lambda_min;
        receiver_2=max(singular./(singular.^2+variance));
        backward=norm(system*receiver-Hs','fro')/(lambda_max*receiver_f+norm(Hs,'fro'));
        reference_error=NaN;mse_reference_error=NaN;condition_error=NaN;
        if reference_check
            reference=right*((singular./(singular.^2+variance)).*left');
            reference_error=norm(receiver-reference,'fro')/max(norm(reference,'fro'),realmin);
            reference_mse=(norm(reference*H-identity,'fro')^2+variance*norm(reference,'fro')^2)/dimension;
            mse_reference_error=abs(direct-reference_mse);
            measured_condition=cond(system,2);condition_error=abs(condition_A-measured_condition)/measured_condition;
        end
        result.Values(point,:)=[direct,consistent,penalty,reg,omit,cross,noise,full_mse,condition_A,lambda_min,lambda_max, ...
            receiver_2,receiver_f,amplified_2,amplified_f,gain_f,gain_2,backward,reconstruction,penalty_error, ...
            reference_error,mse_reference_error,condition_error];
    end
end