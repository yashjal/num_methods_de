function [y_alias, y_anti] = fourier_kdv_rhs(phi_hat, L)
% computes the fourier coefiicients of rhs of KdV
% given fourier coefficients of phi (initial condition function)
    n = size(phi_hat,2);
    x = linspace(-L/2,L/2,n+1); x(end) = [];
    if mod(n,2) == 0    % even number of modes
        ik_cubed = 1i * ((2*pi/L) * [0:n/2-1 0 -n/2+1:-1]).^3;
        ik = 1i * (2*pi/L) * [0:n/2-1 0 -n/2+1:-1];
    else
        ik_cubed = 1i * ((2*pi/L) * [0:(n-1)/2 -(n-1)/2:-1]).^3;
        ik = 1i .* (2*pi/L) .* [0:(n-1)/2 -(n-1)/2:-1];
    end
    y_alias = ik_cubed .* phi_hat - (3*ik) .* fft( ifft(phi_hat).^2 );
    y_anti = ik_cubed .* phi_hat - (3*ik) .* fft( ifft(anti_alias(phi_hat,phi_hat)) );
    
    figure(20); clf;
    plot(x,ifft(phi_hat).^2,'b.-'); hold on;
    plot(x,ifft(anti_alias(phi_hat,phi_hat)),'r.-');
    xlabel('x'); ylabel('(\phi(x))^{2}');
    title(['Approximating \phi^{2} with and without anti-aliasing for n=',num2str(n)]);
    legend('aliased approximate','anti-aliased approximate');
    pause
    %y_anti = 0;
end