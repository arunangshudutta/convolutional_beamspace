function filter_responce_plot(h_coeff)
    % plot of frequency responce of filter
    % h_coeff: filter coefficients
    
    [h_rps_plt,w_rng_plt] = freqz(h_coeff,1,512);
    figure
    plot(w_rng_plt/pi, 20*log10(abs(h_rps_plt)) ); grid on;
    xlabel('Normalized Frequency (\omega/\pi)')
    ylabel('Magnitude (dB)')
    title('Responce of filter')
end