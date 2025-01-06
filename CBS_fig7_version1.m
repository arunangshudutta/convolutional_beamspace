clc; clear;
format long 

N_sim = 5000; % number of montecarlo sim

SNR_dB_vec = -6:2:6;

doa = [-0.5, 0.5]*pi/180; % source doa in radians
D = length(doa); % number of sources
N = 99; % number of sensors
num_snapshot = 100; % number of samples/snapshots 

M = 12; % decimation ratio

% low pass FIR filter parameters
L = 16; % filter length (order of FIR filter: L-1)
pass_ed = pi/(2*M); % passband edge in radians/second
stop_ed = (3*pi)/(2*M); % stopband edge in radians/second

freq_points = [0, pass_ed, stop_ed, pi]; % Normalized frequency points in radians/second
freq_points = freq_points/pi; % (w/pi) input to matlab func firpm, range [0,1]

amp_points = [1, 1, 0, 0]; % amplitude at frequency points

% filter coefficients
h_coeff = firpm(L-1, freq_points, amp_points); 

% plot of frequency responce of filter
%filter_responce_plot(h_coeff);

%----Toplitz matrix----%
H = toplitz_matrix(N,L,h_coeff);


%----Simulation start----%

err_1poly = zeros(length(SNR_dB_vec),1); err_noncoh = err_1poly; err_coh = err_1poly;
for k_snr = 1:length(SNR_dB_vec)
    SNR_dB = SNR_dB_vec(k_snr);
    
    disp(SNR_dB);

    err_1poly_temp = zeros(N_sim,1); err_noncoh_temp = err_1poly_temp; 
    err_coh_temp = err_1poly_temp;

    for k_sim=1:N_sim
    
        %----Received signal----%
        x = sensor_array_output(D,doa,N,SNR_dB,num_snapshot);
        
        % steady state output samples
        y = H*x;
        
        % covariance calculated using element space
        % R_ele = x*(x')/num_snapshot;
        
        % doa estimates and RMSE usig coherent method
        doa_est_coh = doa_estimate_coherent(y,M,D);
        err_coh_temp(k_sim) = sum((doa.' - doa_est_coh).^2);

        % doa estimates and RMSE usig noncoherent method
        doa_est_noncoh = doa_estimate_noncoherent(y,M,D);
        err_noncoh_temp(k_sim) = sum((doa.' - doa_est_noncoh).^2);

        % doa estimates and RMSE usig 1 polyphase component
        doa_est_1poly = doa_estimate_1poly(y,M,D);
        err_1poly_temp(k_sim) = sum((doa.' - doa_est_1poly).^2);
        

    end
    err_coh(k_snr) = sqrt(sum(err_coh_temp)/(N_sim*D));
    err_noncoh(k_snr) = sqrt(sum(err_noncoh_temp)/(N_sim*D));
    err_1poly(k_snr) = sqrt(sum(err_1poly_temp)/(N_sim*D));
end

%% 

figure
semilogy(SNR_dB_vec,err_coh,'-*',LineWidth=1.5)
hold on
semilogy(SNR_dB_vec,err_noncoh,'-^',LineWidth=1.5)
hold on
semilogy(SNR_dB_vec,err_1poly,'-o',LineWidth=2)
grid on; 
legend('Coherent', 'Noncoherent', '1 polyphase component')
ylabel('RMSE');xlabel('SNR in dB')

%% Functions


function R_est = est_covarince_trun(y,N_bs)
    num_snapshot = size(y,2);
    v_vec_temp = y(1:N_bs,:); % truncated vector
    R_est = v_vec_temp*(v_vec_temp')/num_snapshot;
end


function doa_est_coh = doa_estimate_coherent(y,M,D)

    J = ceil(size(y,1)/M); %ceil((N-L+1)/M);
    num_snapshot = size(y,2);
    R_dec_avg = zeros(J,J);
    for k=1:M
        v_vec_temp = y(k:M:end,:); % decimated vector
        % estimated covariance matrix
        R_dec_temp = (v_vec_temp*(v_vec_temp'))/num_snapshot;
        R_dec_avg = R_dec_avg + R_dec_temp;
    end
    R_dec_avg = R_dec_avg/M; % averaged estimate

    w_est_dec = sort(rootmusic(R_dec_avg,D)/M); 
    doa_est_coh = asin(w_est_dec/pi);

end

function doa_est_noncoh = doa_estimate_noncoherent(y,M,D)

    num_snapshot = size(y,2);
    doa_est_noncoh = zeros(D,1);
    for k=1:M
        v_vec_temp = y(k:M:end,:); % decimated vector
        % estimated covariance matrix
        R_dec_temp = (v_vec_temp*(v_vec_temp'))/num_snapshot;
        % estimating doa
        w_est_dec_temp = sort(rootmusic(R_dec_temp,D)/M);
        doa_est_noncoh = doa_est_noncoh + asin(w_est_dec_temp/pi);
    end
    doa_est_noncoh = doa_est_noncoh/M;
end

function doa_est_1poly = doa_estimate_1poly(y,M,D)

    num_snapshot = size(y,2);

    v_vec_temp = y(1:M:end,:); % decimated vector
    R_dec_temp = (v_vec_temp*(v_vec_temp'))/num_snapshot;
    % estimating doa
    w_est_dec_temp = sort(rootmusic(R_dec_temp,D)/M);
    doa_est_1poly = asin(w_est_dec_temp/pi);

end

function H = toplitz_matrix(N,L,h_coeff)
    H = zeros(N-L+1, N);
    h_coeff_flip = flip(h_coeff);
    for k=1:(N-L+1)
        h = [zeros(1,(k-1)),h_coeff_flip,zeros(1,(N-L- (k-1)))];
        H(k,:) = h;
    end
end

function x = sensor_array_output(D,doa,N,SNR_dB,num_snapshot)
    % source signals
    c = randn(D,num_snapshot); % E(c) = 0
    % Vandermonde matrix of steering vectors
    A = zeros(N,D);
    N_vec = 0:(N-1); N_vec = N_vec.';
    for k=1:D
        A(:,k) = exp(1i*pi*sin(doa(k))*N_vec);
    end
    x_nless = A*c; % noiseless received signal
    e = awgn(x_nless,SNR_dB); % noise 
    x = x_nless + e; % sensor array output
end


function filter_responce_plot(h_coeff)
    % plot of frequency responce of filter
    [h_rps_plt,w_rng_plt] = freqz(h_coeff,1,512);
    figure
    plot(w_rng_plt/pi, 20*log10(abs(h_rps_plt)) ); grid on;
    xlabel('Normalized Frequency (\omega/\pi)')
    ylabel('Magnitude (dB)')
    title('Responce of filter used')
end