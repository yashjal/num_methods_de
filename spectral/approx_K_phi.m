% 1.2 hw2 pseudo-spectral approx of K[phi(.)]

L = 60;
c = 1;
err1 = []; err2 = []; errinf = [];
x_plot = 4:2:500;
n_tol = []; L_plot = []; c_plot = []; n_peaks = [];
err_tol = 1e-9;

%for L=40:10:100
%for c=1.0:0.25:3.0
for n=4:2:500
    N = n+1000;

    x = linspace(-L/2,L/2,n+1); x(end) = [];
    x_fine = linspace(-L/2,L/2,N+1); x_fine(end) = [];

    f = intial_cond_deriv_soliton(x,c);
    f_fine = intial_cond_deriv_soliton(x_fine,c);

    phi_hat = fft(intial_cond_soliton(x,c));
    [K_hat,~] = fourier_kdv_rhs(phi_hat,L);
    K_hat = K_hat ./ n;

    if mod(n,2) == 0
        K_approx= N * ifft([K_hat(1:n/2-1) -0.5*K_hat(n/2) zeros(1,N-n-1) 0.5*K_hat(n/2) K_hat(n/2+1:n)]);
    else
        K_approx= N * ifft([K_hat(1:(n-1)/2) zeros(1,N-n) K_hat((n+1)/2:n)]);
    end

    err1 = [err1, e1_rel_err(f_fine,K_approx,L)];
    err2 = [err2, e2_rel_err(f_fine,K_approx,L)];
    errinf = [errinf, inf_rel_err(f_fine,K_approx)];

%     if (e1_rel_err(f_fine,K_approx,L) < err_tol || inf_rel_err(f_fine,K_approx) < err_tol || e2_rel_err(f_fine,K_approx,L) <err_tol)
%         n_tol = [n_tol, n];
%         %L_plot = [L_plot,L];
%         c_plot = [c_plot, c];
%         
%         max_f = max(f_fine);
%         half_f = 0.5*max_f;
%         [min_diff,index_half_f] = min(abs(half_f-f));
%         n_peak = n-(2*(n-index_half_f));
%         n_peaks = [n_peaks, n_peak];
%         break;
%     end

    if n == 50
        figure(2); clf;
        plot(x_fine, f_fine,'r.-'); hold on;
        plot(x_fine, K_approx, 'b.-');
        title(['Approximating K[\phi_{sol}(.,t=0)] using Fourier series of \phi_{sol}(x,t=0) for n=',num2str(n)]);
        xlabel('x'); ylabel('y');
        legend('K[\phi_{sol}(x,t=0)] actual function','\phi^{*}(x) approximate');

        figure(3); clf;
        plot(x_fine, abs(f_fine-K_approx),'k.-');
        title(['Absolute error in approximating K[\phi_{sol}(.,t=0)] using Fourier series of \phi_{sol}(x,t=0) for n=',num2str(n)]);
        xlabel('x'); ylabel('y');

        %figure(4); clf;
        %plot(x, abs(K_hat),'k.-');
    end
end
%end

figure(1); clf;
semilogy(x_plot, err1, 'r.-'); hold on;
semilogy(x_plot, err2, 'b.-');
semilogy(x_plot, errinf, 'k.-');
title('Error in approximating K[\phi_{sol}(.,t=0)]=d_{t}\phi_{sol}(x,t=0) using truncated Fourier series of \phi_{sol}(x,t=0)');
xlabel('n number of Fourier modes'); ylabel('\epsilon relative error');
legend('\epsilon_{1}','\epsilon_{2}','\epsilon_{\infty}');

% figure(5); clf;
% plot(L_plot, n_tol, 'ko-'); 
% xlabel('L length of domain interval'); ylabel('N number of modes to achieve 1e-9 accuracy');

% figure(6); clf;
% plot(c_plot, n_tol, 'ko-'); 
% xlabel('c'); ylabel('N number of modes to achieve 1e-9 accuracy');