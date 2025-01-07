function x = sensor_array_output(doa,N,SNR_dB,num_snapshot)
    
    % doa: DOAs in radians
    % N: # of sensors
    % SNR_dB: SNR in dB
    % num_snapshots: # of snapshots or samples of received vector

    D = length(doa); % number of sources;

    % source signals
    c = randn(D,num_snapshot); % E(c) = 0
    % power of source signal = 1

    % Vandermonde matrix of steering vectors
    A = zeros(N,D);
    N_vec = 0:(N-1); N_vec = N_vec.';
    for k=1:D
        A(:,k) = exp(1i*pi*sin(doa(k))*N_vec);
    end

    x_nless = A*c; % noiseless received signal
    % Power of received signal = D*(power of source signal) = D
    
    signal_power = D;
    noise_power = signal_power/(10^(0.1*SNR_dB));

    % AWGN noise
    e = sqrt(noise_power/2)*randn(N,num_snapshot) + 1i*sqrt(noise_power/2)*randn(N,num_snapshot);

    x = x_nless + e; % sensor array output
end