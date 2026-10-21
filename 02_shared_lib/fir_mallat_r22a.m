function [P,dP,meta]=fir_mallat_r22a(dimension,alpha)
    assert(dimension>=3&&mod(dimension,1)==0,'R22A:Dimension','Integer Abid axis required.');
    rotation=[cos(alpha),-sin(alpha);sin(alpha),cos(alpha)];
    derivative=[-sin(alpha),-cos(alpha);cos(alpha),-sin(alpha)];
    Q=[1,1;1,-1]/sqrt(2);
    E=zeros(2,2,2);dE=E;
    for phase=1:2
        vector=rotation(:,phase);tangent=derivative(:,phase);
        E(:,:,phase)=vector*vector'*Q;
        dE(:,:,phase)=(tangent*vector'+vector*tangent')*Q;
    end
    h=reshape(E(:,1,:),1,[]);g=reshape(E(:,2,:),1,[]);
    dh=reshape(dE(:,1,:),1,[]);dg=reshape(dE(:,2,:),1,[]);
    [P,dP,levels]=cascade_r22a(dimension,h,g,dh,dg);
    degree_residual=norm(E(:,:,1)*E(:,:,2)','fro');
    paraunitary=norm(E(:,:,1)*E(:,:,1)'+E(:,:,2)*E(:,:,2)'-eye(2),'fro')+2*degree_residual;
    meta=struct('h',h,'g',g,'dh',dh,'dg',dg,'E',E,'dE',dE,'Levels',levels, ...
        'Depth',sum(mod(dimension./2.^(0:floor(log2(dimension))),2)==0), ...
        'ParaunitarityResidual',paraunitary,'DCResidual',abs(sum(h)-sqrt(2)), ...
        'HighpassDCResidual',abs(sum(g)),'CQFResidual',norm(g-(-1).^(0:3).*fliplr(h)), ...
        'OrthogonalityResidual',norm(P'*P-eye(dimension),'fro'), ...
        'FilterLength',4,'Alpha',alpha,'IntrinsicDimension',1);
end

function [P,dP,levels]=cascade_r22a(dimension,h,g,dh,dg)
    if mod(dimension,2)==1
        P=eye(dimension);dP=zeros(dimension);levels=zeros(1,dimension);return
    end
    [coarse,dcoarse,previous_levels]=cascade_r22a(dimension/2,h,g,dh,dg);
    low=zeros(dimension,dimension/2);high=low;dlow=low;dhigh=low;
    for column=1:dimension/2
        pulse=zeros(1,dimension);pulse(2*column-1)=1;
        filters={h,g,dh,dg};outputs=cell(1,4);
        for filter_index=1:4
            signal=conv(pulse,filters{filter_index});
            outputs{filter_index}=accumarray(mod((0:numel(signal)-1).',dimension)+1,signal(:),[dimension,1]);
        end
        low(:,column)=outputs{1};high(:,column)=outputs{2};
        dlow(:,column)=outputs{3};dhigh(:,column)=outputs{4};
    end
    P=[low*coarse,high];dP=[dlow*coarse+low*dcoarse,dhigh];
    levels=[previous_levels+1,ones(1,dimension/2)];
end