

function music_pseudo_spectrum(D,N,R)

    [eig_vec_mat,temp]=eig(R); %Find the eigenvalues and eigenvectors of R
    
    % estimated noise subspace
    noise_subspace = eig_vec_mat(:,1:(N-D)) ; % N-D coulmns (associated with smallest N-D eigenvalues)
    
    % DOA sweep search range -90 to 90 degrees
    theta_sweep = -90:0.5:90;
    theta_sweep_r = theta_sweep*pi/180;
    
    % calculating pseudospectrum
    N_vec = 0:(N-1); N_vec = N_vec.';
    music_spec = zeros(1,length(theta_sweep));
    for k = 1:length(theta_sweep)
        theta_temp = theta_sweep_r(k);
        a = exp(1i*pi*sin(theta_temp)*N_vec);
        music_spec(k) = 1/abs((a')*noise_subspace*(noise_subspace')*a);
    end
    
    music_spec_dB=10*log10(music_spec/max(music_spec)); 
    
    figure
    plot(theta_sweep,music_spec_dB,'-k')
    xlabel('angle \theta (radians)')
    ylabel('pseudospectrum P(\theta) /dB')
    title('DOA estimation based on MUSIC algorithm (normalized + dB scale)')
    grid on

end