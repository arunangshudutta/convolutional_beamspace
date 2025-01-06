clc; clear;
format long 

N_sim = 1000; % number of montecarlo sim
N_bs_vec = 7:7:84; 
SNR_dB = 0;

doa = [-5, 5]*pi/180; % source doa in radians
D = length(doa); % number of sources
N = 99; % number of sensors
num_snapshot = 100; % number of samples/snapshots 

M = 4; % decimation ratio

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

err_ele = zeros(length(N_bs_vec),1); err_trun = err_ele; err_dec = err_ele;
for k_bs = 1:length(N_bs_vec)
    N_bs = N_bs_vec(k_bs);
    disp(N_bs);

    err_ele_temp = zeros(N_sim,1); err_trun_temp = err_ele_temp; 
    err_dec_temp = err_ele_temp;
    for k_sim=1:N_sim
    
        %----Received signal----%
        x = sensor_array_output(D,doa,N,SNR_dB,num_snapshot);
        
        % steady state output samples
        y = H*x;
        
        % covariance calculated using element space
        R_ele = x*(x')/num_snapshot;
        
        % covariance calculated using decimated vector
        R_dec_avg = avaraged_covarince_dec(y,M);
        
        % covariance calculated using truncated vector
        R_trun_avg = est_covarince_trun(y,N_bs);
        
        % Root-MUSIC and RMSE % need to change this part to acount for unknown #D 
        w_est_ele = sort(rootmusic(R_ele,D)); 
        doa_est_ele = asin(w_est_ele/pi); % doa estimates using element space
        err_ele_temp(k_sim) = sum((doa.' - doa_est_ele).^2);
        
        w_est_dec = sort(rootmusic(R_dec_avg,D)/M); 
        doa_est_dec = asin(w_est_dec/pi); % doa estimates using CBS decimated vector
        err_dec_temp(k_sim) = sum((doa.' - doa_est_dec).^2);
        
        w_est_trun = sort(rootmusic(R_trun_avg,D)); 
        doa_est_trun = asin(w_est_trun/pi); % doa estimates using CBS decimated vector
        err_trun_temp(k_sim) = sum((doa.' - doa_est_trun).^2);
    end
    err_ele(k_bs) = sqrt(sum(err_ele_temp)/(N_sim*D));
    err_dec(k_bs) = sqrt(sum(err_dec_temp)/(N_sim*D));
    err_trun(k_bs) = sqrt(sum(err_trun_temp)/(N_sim*D));
end

%% 

figure
semilogy(N_bs_vec,err_ele,LineWidth=2)
hold on
semilogy(N_bs_vec,err_dec,'-^',LineWidth=2)
hold on
semilogy(N_bs_vec,err_trun,'-o',LineWidth=2)
grid on; xlim([0 90]);
legend('element space', 'CBS, decimated M=4', 'CSB, truncated')
ylabel('RMSE');xlabel('Beamspace size N_{bs}')
title(['RMSE of truncated CBS, decimated CBS, and element-space; ' 'SNR = ',num2str(SNR_dB),'dB; '])
%% Functions


function R_est = est_covarince_trun(y,N_bs)
    num_snapshot = size(y,2);
    v_vec_temp = y(1:N_bs,:); % truncated vector
    R_est = v_vec_temp*(v_vec_temp')/num_snapshot;
end


function R_dec_avg = avaraged_covarince_dec(y,M)
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