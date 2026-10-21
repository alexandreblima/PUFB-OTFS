function result=r22dcb_math(action,varargin)
    switch action
        case 'receiver'
            H=varargin{1};variance=varargin{2};factor=chol(H'*H+variance*eye(size(H,2)));
            result=factor\(factor'\H');
        case 'bootstrap'
            result=bootstrap_delta(varargin{:});
        case 'crossing'
            result=crossing(varargin{:});
        case 'wilson'
            result=wilson(varargin{:});
        case 'read'
            result=readtable(varargin{1},'Delimiter',',','VariableNamingRule','preserve','TextType','string');
        case 'arrayhash'
            array=varargin{1};if isa(array,'uint8'),bytes=array(:);else,bytes=typecast(array(:),'uint8');end
            digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
            result=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
        case 'channelmean'
            matrix=varargin{1};points=varargin{2};
            if isempty(points),result=(mean(matrix(:,1:3),2)+mean(matrix(:,4:6),2)+mean(matrix(:,7:8),2))/3;
            else,result=mean(matrix(:,points),2);end
        case 'selftest'
            result=selftest();
        otherwise
            error('R22DCB:Action','Unknown math action.');
    end
end

function result=bootstrap_delta(values,channels_per_velocity,velocity_count,seed,replicates)
    values=reshape(values,channels_per_velocity,velocity_count);stream=RandStream('mt19937ar','Seed',seed);
    means=zeros(replicates,velocity_count);
    for velocity=1:velocity_count
        draws=randi(stream,channels_per_velocity,[channels_per_velocity,replicates]);current=values(:,velocity);
        means(:,velocity)=mean(current(draws),1).';
    end
    samples=sort(mean(means,2));
    result=struct('Mean',mean(values,'all'),'Median',median(values(:)),'SD',std(values(:)), ...
        'Standardized',mean(values,'all')/max(std(values(:)),eps), ...
        'Lower',samples(ceil(.025*replicates)),'Upper',samples(ceil(.975*replicates)));
end

function result=crossing(snr,rate,target,eligible)
    result=struct('SNR',NaN,'LeftIndex',NaN,'RightIndex',NaN,'Status','TARGET_NOT_RESOLVED');
    for point=1:numel(snr)-1
        if eligible(point)&&eligible(point+1)&&isfinite(rate(point))&&isfinite(rate(point+1)) ...
                &&rate(point)>0&&rate(point+1)>0&&rate(point)>=target&&rate(point+1)<=target&&rate(point)>rate(point+1)
            weight=(log10(target)-log10(rate(point)))/(log10(rate(point+1))-log10(rate(point)));
            result.SNR=snr(point)+weight*(snr(point+1)-snr(point));result.LeftIndex=point;result.RightIndex=point+1;
            result.Status='RESOLVED';return
        end
    end
end

function interval=wilson(errors,bits,alpha)
    if nargin<3,alpha=.05;end
    rate=errors/bits;critical=sqrt(2)*erfcinv(alpha);denominator=1+critical^2/bits;
    center=(rate+critical^2/(2*bits))/denominator;
    radius=critical*sqrt(rate*(1-rate)/bits+critical^2/(4*bits^2))/denominator;
    interval=[max(0,center-radius),min(1,center+radius)];
end

function result=selftest()
    stream=RandStream('mt19937ar','Seed',131001);H=randn(stream,12)+1i*randn(stream,12);variance=.1;
    receiver=r22dcb_math('receiver',H,variance);direct=(H'*H+variance*eye(12))\H';
    assert(norm(receiver-direct,'fro')<1e-11);
    resolved=crossing(0:3,[.1,.01,.001,.0001],.003,true(1,4));
    assert(abs(resolved.SNR-(1+log10(.01/.003)))<1e-12);
    unresolved=crossing(0:3,[.1,.02,.01,.005],.001,true(1,4));assert(isnan(unresolved.SNR));
    zero=crossing(0:2,[.01,.002,0],.001,true(1,3));assert(isnan(zero.SNR));
    interval=wilson(0,10000);assert(interval(1)==0&&interval(2)>0);
    result=struct('PASS',true,'ReceiverSolveError',norm(receiver-direct,'fro'),'NoExtrapolationPASS',true);
end