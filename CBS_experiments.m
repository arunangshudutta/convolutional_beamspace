clc; clear;

%% Understanding CBS using MUSIC-algo pseudo-spectrum

w = 0.5;
doa = asin(w/pi);
doa_deg = doa*180/pi;
disp(doa_deg)

%doa = [5]*pi/180; % source doa in radians
D = length(doa); % number of sources
N = 99; % number of sensors
num_snapshot = 100; % number of samples/snapshots 
SNR_dB = 0;


% low pass FIR filter parameters
L = 16; % filter length (order of FIR filter: L-1)
M = 4;
pass_ed = pi/(2*M); % passband edge in radians/second
stop_ed = (3*pi)/(2*M); % stopband edge in radians/second

freq_points = [0, pass_ed, stop_ed, pi]; % Normalized frequency points in radians/second
freq_points = freq_points/pi; % (w/pi) input to matlab func firpm, range [0,1]

amp_points = [1, 1, 0, 0]; % amplitude at frequency points

% filter coefficients
h_coeff = firpm(L-1, freq_points, amp_points); 

% plot of frequency responce of filter
% filter_responce_plot(h_coeff);

%----Toplitz matrix----%
H = filter_toplitz_matrix(N,L,h_coeff);

x = sensor_array_output(doa,N,SNR_dB,num_snapshot);

y = H*x;

R_ele = x*(x')/num_snapshot;

R_cbs = y*(y')/num_snapshot;

R = R_cbs; 

%%

[eig_vec_mat,temp]=eig(R); %Find the eigenvalues and eigenvectors of R
N = size(R,1);
% estimated noise subspace
noise_subspace = eig_vec_mat(:,1:(N-D)) ; % N-D coulmns (associated with smallest N-D eigenvalues)



[h_rps_plt,theta_sweep_r] = freqz(h_coeff,1,1024*10);


% calculating pseudospectrum
N_vec = 0:(N-1); N_vec = N_vec.';
music_spec = zeros(1,length(theta_sweep_r));
for k = 1:length(theta_sweep_r)
    theta_temp = theta_sweep_r(k);
    a = exp(1i*pi*sin(theta_temp)*N_vec);
    music_spec(k) = 1/abs((a')*noise_subspace*(noise_subspace')*a);
end

music_spec_dB=10*log10(music_spec/max(music_spec)); 

figure
plot(theta_sweep_r, 20*log10(abs(h_rps_plt)));
hold on
plot(pi*sin(theta_sweep_r),music_spec_dB,'-k')
hold on
plot(pass_ed,0,'*', LineWidth=2)
hold on
plot(stop_ed,0,'*', LineWidth=2)
legend('filter', 'doa spectrum (w = pi*sin(\theta))', 'pass band edge', 'stop band edge')
xlabel('w radians')
ylabel('Magnitude (dB)')
grid on

% figure
% plot(theta_sweep_r*180/pi, 20*log10(abs(h_rps_plt)));
% hold on
% plot(180*sin(theta_sweep_r),music_spec_dB,'-k')
% hold on
% plot(pass_ed*180/pi,0,'*', LineWidth=2)
% hold on
% plot(stop_ed*180/pi,0,'*', LineWidth=2)
% legend('filter', 'doa spectrum', 'pass band edge', 'stop band edge')
% xlabel('degrees')
% ylabel('Magnitude (dB)')
% grid on


%music_pseudo_spectrum(D,(N-L+1),R_cbs);