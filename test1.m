clc;clear;

% signal power
w= pi/6;
N=10000;

x = exp(1i*w*[1:N]);

pow_sig = (x*(x'))/N;

% noise power

N=10000000;
P_n_dB = 34; %dB
P_n = 10^(0.1*P_n_dB);
n = sqrt(P_n/2)*randn(1,N) + 1i*sqrt(P_n/2)*randn(1,N);

pow_noise = (n*(n'))/N;

pow_noise_dB = 10*log10(pow_noise);


%%

clc; clear;
doa = [-5, 5, 10]*pi/180;
D = 3; num_snapshot = 10000;
SNR_dB = 0;
N = 99;


% source signals
c = randn(D,num_snapshot); % E(c) = 0

% source signal power 
pow_c = (c(1,:)*(c(1,:)'))/num_snapshot; % varince of c = 1

% Vandermonde matrix of steering vectors
A = zeros(N,D);
N_vec = 0:(N-1); N_vec = N_vec.';
for k=1:D
    A(:,k) = exp(1i*pi*sin(doa(k))*N_vec);
end
x_nless = A*c; % noiseless received signal

% received signal power
k = 50;
pow_xnless = (x_nless(k,:)*(x_nless(k,:)'))/num_snapshot; % power = D*(power of c) = 2*1

signal_power = D; 
noise_power = signal_power/(10^(0.1*SNR_dB));

e = sqrt(noise_power/2)*randn(N,num_snapshot) + 1i*sqrt(noise_power/2)*randn(N,num_snapshot);

%e = awgn(x_nless,SNR_dB); % noise 


% noise power
pow_e = (e(k,:)*(e(k,:)'))/num_snapshot;

x = x_nless + e; % sensor array output

