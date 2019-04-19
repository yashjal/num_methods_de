% 1.1 hw2 initial condition approximation error

L = 60;
c = 1;
err1_sol = []; err2_sol = []; errinf_sol = [];
err1_tr = []; err2_tr = []; errinf_tr = [];
err1_saw = [];err2_saw = []; errinf_saw = [];
x_plot = 4:2:250;

for n=4:2:250
    N = n+1000;

    x = linspace(-L/2,L/2,n+1); x(end) = [];
    x_fine = linspace(-L/2,L/2,N+1); x_fine(end) = [];

    f_sol = intial_cond_soliton(x,c);
    f_tr = intial_cond_tr(x);
    f_saw = intial_cond_sawtooth(x);

    f_sol_fine = intial_cond_soliton(x_fine,c);
    f_tr_fine = intial_cond_tr(x_fine);
    f_saw_fine = intial_cond_sawtooth(x_fine);

    if mod(n,2) == 0
        fhat_sol = fft(f_sol)./n;
        f_sol_approx = N * ifft([fhat_sol(1:n/2-1) 0.5*fhat_sol(n/2) zeros(1,N-n-1) 0.5*fhat_sol(n/2) fhat_sol(n/2+1:n)]);
        fhat_tr = fft(f_tr)./n;
        f_tr_approx = N * ifft([fhat_tr(1:n/2-1) 0.5*fhat_tr(n/2) zeros(1,N-n-1) 0.5*fhat_tr(n/2) fhat_tr(n/2+1:n)]);
        fhat_saw = fft(f_saw)./n;
        f_saw_approx = N * ifft([fhat_saw(1:n/2-1) 0.5*fhat_saw(n/2) zeros(1,N-n-1) 0.5*fhat_saw(n/2) fhat_saw(n/2+1:n)]);        
    else
        f_sol_approx = interpft(f_sol,N);
        f_tr_approx = interpft(f_tr,N);
        f_saw_approx = interpft(f_saw,N);
    end
    
    err1_sol = [err1_sol, e1_rel_err(f_sol_fine,f_sol_approx,L)];
    err2_sol = [err2_sol, e2_rel_err(f_sol_fine,f_sol_approx,L)];
    errinf_sol = [errinf_sol, inf_rel_err(f_sol_fine,f_sol_approx)];
    
    err1_tr = [err1_tr, e1_rel_err(f_tr_fine,f_tr_approx,L)];
    err2_tr = [err2_tr, e2_rel_err(f_tr_fine,f_tr_approx,L)];
    errinf_tr = [errinf_tr, inf_rel_err(f_tr_fine,f_tr_approx)];
    
    err1_saw = [err1_saw, e1_rel_err(f_saw_fine,f_saw_approx,L)];
    err2_saw = [err2_saw, e2_rel_err(f_saw_fine,f_saw_approx,L)];
    errinf_saw = [errinf_saw, inf_rel_err(f_saw_fine,f_saw_approx)];
    
    if n == 100
        figure(4); clf;
        plot(x_fine, f_sol_fine,'b.-'); hold on;
        plot(x_fine, f_sol_approx,'r.-');
        title(['Approximating the soliton wave using truncated Fourier series with n=',num2str(n)]);
        xlabel('x'); ylabel('y');
        legend('\phi(x,t=0) actual function','\phi^{*}(x) approximation');
        
        figure(5); clf;
        plot(x_fine, abs(f_tr_fine-f_tr_approx),'b.-'); hold on;
        %plot(x_fine, f_tr_approx,'r.-');
        %title(['Approximating the triangle wave f(x) using truncated Fourier series with n=',num2str(n)]);
        title(['Error in approximating the triangle wave f(x) using truncated Fourier series with n=',num2str(n)]);
        xlabel('x'); ylabel('y');
        %legend('f(x) actual function','\phi(x) approximation');
        
        figure(6); clf;
        plot(x_fine, abs(f_saw_fine-f_saw_approx),'b.-'); hold on;
        %plot(x_fine, f_saw_approx,'r.-');
        %title(['Approximating the sawtooth wave g(x) using truncated Fourier series with n=',num2str(n)]);
        title(['Error in approximating the sawtooth wave g(x) using truncated Fourier series with n=',num2str(n)]);
        xlabel('x'); ylabel('y');
        %legend('g(x) actual function','\phi(x) approximation');
    end
end

figure(1); clf;
semilogy(x_plot, err1_sol, 'r.-'); hold on;
semilogy(x_plot, err2_sol, 'b.-');
semilogy(x_plot, errinf_sol, 'k.-');
title('Error in approximating the soliton wave initial condition \phi(x,t=0)=1/2 sech^{2}(1/2 x) using truncated Fourier series');
xlabel('n number of Fourier modes'); ylabel('\epsilon relative error');
legend('\epsilon_{1}','\epsilon_{2}','\epsilon_{\infty}');

figure(2); clf;
semilogy(x_plot, err1_tr, 'r.-'); hold on;
semilogy(x_plot, err2_tr, 'b.-');
semilogy(x_plot, errinf_tr, 'k.-');
title('Error in approximating the triangle wave initial condition f(x) using truncated Fourier series');
xlabel('n number of Fourier modes'); ylabel('\epsilon relative error');
legend('\epsilon_{1}','\epsilon_{2}','\epsilon_{\infty}');

figure(3); clf;
plot(x_plot, err1_saw, 'r.-'); hold on;
plot(x_plot, err2_saw, 'b.-');
plot(x_plot, errinf_saw, 'k.-');
title('Error in approximating the sawtooth wave initial condition g(x) using truncated Fourier series');
xlabel('n number of Fourier modes'); ylabel('\epsilon relative error');
legend('\epsilon_{1}','\epsilon_{2}','\epsilon_{\infty}');