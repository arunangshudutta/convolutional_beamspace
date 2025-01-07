function H = filter_toplitz_matrix(N,L,h_coeff)
    % Toplitz matrix: multiplied with the sensor input to get steady state responce of filter

    % N: # of sensors
    % L: filter length
    % h_coeff: FIR filter coefficients

    H = zeros(N-L+1, N);
    h_coeff_flip = flip(h_coeff);
    for k=1:(N-L+1)
        h = [zeros(1,(k-1)),h_coeff_flip,zeros(1,(N-L- (k-1)))];
        H(k,:) = h;
    end
    
end