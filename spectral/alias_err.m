% 1.3 hw2 aliasing error

L = 60;
c = 1;
n = 20;
N = n + 1000;

x = linspace(-L/2,L/2,n+1); x(end) = [];
x_fine = linspace(-L/2,L/2,N+1); x_fine(end) = [];
f_fine = intial_cond_deriv_soliton(x_fine,c);

phi_hat = fft(intial_cond_soliton(x,c));
[K_hat_alias, K_hat_anti] = fourier_kdv_rhs(phi_hat,L);
K_hat_alias = K_hat_alias ./ n;
K_hat_anti =  K_hat_anti ./ n;

f_hat = fft(intial_cond_tr(x));
g_hat = fft(intial_cond_sawtooth(x));

[f_hat_alias, f_hat_anti] = fourier_kdv_rhs(f_hat,L);
f_hat_alias = f_hat_alias ./ n;
f_hat_anti =  f_hat_anti ./ n;
[g_hat_alias, g_hat_anti] = fourier_kdv_rhs(g_hat,L);
g_hat_alias = g_hat_alias ./ n;
g_hat_anti =  g_hat_anti ./ n;

if mod(n,2) == 0
    K_approx_alias = N * ifft([K_hat_alias(1:n/2-1) -0.5*K_hat_alias(n/2) zeros(1,N-n-1) 0.5*K_hat_alias(n/2) K_hat_alias(n/2+1:n)]);
    K_approx_anti = N * ifft([K_hat_anti(1:n/2-1) -0.5*K_hat_anti(n/2) zeros(1,N-n-1) 0.5*K_hat_anti(n/2) K_hat_anti(n/2+1:n)]);
    f_approx_alias = N * ifft([f_hat_alias(1:n/2-1) -0.5*f_hat_alias(n/2) zeros(1,N-n-1) 0.5*f_hat_alias(n/2) f_hat_alias(n/2+1:n)]);
    f_approx_anti = N * ifft([f_hat_anti(1:n/2-1) -0.5*f_hat_anti(n/2) zeros(1,N-n-1) 0.5*f_hat_anti(n/2) f_hat_anti(n/2+1:n)]);
    g_approx_alias = N * ifft([g_hat_alias(1:n/2-1) -0.5*g_hat_alias(n/2) zeros(1,N-n-1) 0.5*g_hat_alias(n/2) g_hat_alias(n/2+1:n)]);
    g_approx_anti = N * ifft([g_hat_anti(1:n/2-1) -0.5*g_hat_anti(n/2) zeros(1,N-n-1) 0.5*g_hat_anti(n/2) g_hat_anti(n/2+1:n)]);
else
    K_approx_alias = N * ifft([K_hat_alias(1:(n-1)/2) zeros(1,N-n) K_hat_alias((n+1)/2:n)]);
    K_approx_anti = N * ifft([K_hat_anti(1:(n-1)/2) zeros(1,N-n) K_hat_anti((n+1)/2:n)]);
    f_approx_alias = N * ifft([f_hat_alias(1:(n-1)/2) zeros(1,N-n) f_hat_alias((n+1)/2:n)]);
    f_approx_anti = N * ifft([f_hat_anti(1:(n-1)/2) zeros(1,N-n) f_hat_anti((n+1)/2:n)]);
    g_approx_alias = N * ifft([g_hat_alias(1:(n-1)/2) zeros(1,N-n) g_hat_alias((n+1)/2:n)]);
    g_approx_anti = N * ifft([g_hat_anti(1:(n-1)/2) zeros(1,N-n) g_hat_anti((n+1)/2:n)]);
end

figure(1); clf;
plot(x_fine,K_approx_alias,'r.-'); hold on;
plot(x_fine,K_approx_anti,'b.-');
xlabel('x'); ylabel('y');
title(['Approximating K[\phi_{sol}(x,t=0)] with and without anti-aliasing for n=',num2str(n)]);
legend('aliased approximate','anti-aliased approximate');

figure(4); clf;
plot(x_fine,abs(K_approx_alias-K_approx_anti),'k.-');
xlabel('x'); ylabel('|K[\phi(x)]-K_{n}[\phi(x)]|');
title(['Absolute error between aliased and anti-aliased approximation of K[\phi_{sol}(x,t=0)] for n=',num2str(n)]);

figure(11); clf;
freq = 2*pi/L * (-n/2:(n/2-1)); % assuming N even
plot(freq,abs(fftshift(K_hat_alias)),'b.-'); hold on;
plot(freq,abs(fftshift(K_hat_anti)),'r.-');
xlabel('freq'); ylabel('$|\hat{K[\phi_{sol}]}|$','Interpreter','latex');
title(['Fourier spectrum of K[\phi_{sol}(x)] for n=',num2str(n)]);
legend('aliased approximate','anti-aliased approximate');

figure(2); clf;
plot(x_fine,f_approx_alias,'r.-'); hold on;
plot(x_fine,f_approx_anti,'b.-');
xlabel('x'); ylabel('y');
title(['Approximating K[f(x)] with and without anti-aliasing for n=',num2str(n)]);
legend('aliased approximate','anti-aliased approximate');

figure(5); clf;
plot(x_fine,abs(f_approx_alias-f_approx_anti),'k.-');
xlabel('x'); ylabel('|K[f(x)]-K_{n}[f(x)]|');
title(['Absolute error between aliased and anti-aliased approximation of K[f(x)] for n=',num2str(n)]);

figure(12); clf;
plot(freq,abs(fftshift(f_hat_alias)),'b.-'); hold on;
plot(freq,abs(fftshift(f_hat_anti)),'r.-');
xlabel('freq'); ylabel('$|\hat{K[f]}|$','Interpreter','latex');
title(['Fourier spectrum of K[f(x)] for n=',num2str(n)]);
legend('aliased approximate','anti-aliased approximate');

figure(3); clf;
plot(x_fine,g_approx_alias,'r.-'); hold on;
plot(x_fine,g_approx_anti,'b.-');
xlabel('x'); ylabel('y');
title(['Approximating K[g(x)] with and without anti-aliasing for n=',num2str(n)]);
legend('aliased approximate','anti-aliased approximate');

figure(6); clf;
plot(x_fine,abs(g_approx_alias-g_approx_anti),'k.-');
xlabel('x'); ylabel('|K[g(x)]-K_{n}[g(x)]|');
title(['Absolute error between aliased and anti-aliased approximation of K[g(x)] for n=',num2str(n)]);

figure(13); clf;
plot(freq,abs(fftshift(g_hat_alias)),'b.-'); hold on;
plot(freq,abs(fftshift(g_hat_anti)),'r.-');
xlabel('freq'); ylabel('$|\hat{K[g]}|$','Interpreter','latex');
title(['Fourier spectrum of K[g(x)] for n=',num2str(n)]);
legend('aliased approximate','anti-aliased approximate');
