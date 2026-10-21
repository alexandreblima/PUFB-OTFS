function [G,info]=generate_channel_r20a(M,N,velocity_mps,seed,protocol)
    L=M*N;table=protocol.eva_table1;delay_ns=table.delays_ns(:).';
    power=10.^(table.relative_powers_db(:).'/10);power=power/sum(power);
    delta_f=protocol.abid_table3.subcarrier_spacing_hz;
    persistent cache_M cache_N cache_delta_f cache_delay_ns cache_Pi
    if isempty(cache_Pi)||isempty(cache_M)||cache_M~=M||cache_N~=N||cache_delta_f~=delta_f ...
            ||~isequal(cache_delay_ns,delay_ns)
        q=(0:L-1).';F=dftmtx(L)/sqrt(L);delay_index=M*(delay_ns*1e-9)*delta_f;
        cache_Pi=complex(zeros(L,L,numel(delay_index)));
        for path=1:numel(delay_index)
            cache_Pi(:,:,path)=F'*diag(exp(-1i*2*pi*q*delay_index(path)/L))*F;
        end
        cache_M=M;cache_N=N;cache_delta_f=delta_f;cache_delay_ns=delay_ns;
    end
    stream=RandStream('mt19937ar','Seed',seed);
    gains=(randn(stream,1,numel(power))+1i*randn(stream,1,numel(power)))./sqrt(2).*sqrt(power);
    theta=-pi+2*pi*rand(stream,1,numel(power));
    carrier=protocol.abid_table3.carrier_frequency_hz;
    T=1/delta_f;
    tau=delay_ns*1e-9;delay_index=M*tau*delta_f;
    doppler_hz=(carrier*velocity_mps/299792458)*cos(theta);
    doppler_index=doppler_hz*N*T;
    n=(0:L-1).';G=complex(zeros(L,L));
    for path=1:numel(power)
        doppler_phase=exp(1i*2*pi*doppler_index(path)*n/L);
        G=G+gains(path)*(cache_Pi(:,:,path).*doppler_phase.');
    end
    info=struct('seed',seed,'velocity_mps',velocity_mps,'path_delays_ns',delay_ns, ...
        'path_powers',power,'path_gains',gains,'theta_rad',theta,'doppler_hz',doppler_hz, ...
        'doppler_index',doppler_index,'delay_index',delay_index, ...
        'frobenius_energy',norm(G,'fro')^2,'dimension',L);
end