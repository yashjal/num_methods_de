N = 225; % # Fourier modes odd
M = N + 1000; % finer grid size
L = 60; % length of domain interval
c = 1; % speed of soliton
x = linspace(-L/2,L/2,N+1)'; x(end) = [];
x_fine = linspace(-L/2,L/2,M+1)'; x_fine(end) = [];
f_fine = intial_cond_soliton(x_fine,c,0); %intial_cond_soliton2(x_fine,c,3*c,0)
T = L/c; % time period
h = 2e-4; % delta t

hlis = 0.0001 * 2.^(0:10);
%err1 = []; err2 = []; err3 = [];
philis = [];

%ylimits=[-0.5 2];
%figure(1); clf
%p = plot(x_fine,f_fine,'linewidth',3);
%xlim([-L/2 L/2])
%ylim(ylimits);
%grid on

% pre-assign values that doesn't change after each time step
k = (2*pi/L) * [0:(N-1)/2 -(N-1)/2:-1]'; % Fourier freq
A = 1i*k.^3;
g = -3i*k;
for h=hlis
invA = ones(N,1)-(2*h/3)*A;
phi_hat_1 = fft(intial_cond_soliton(x,c,h));
phi_hat_0 = fft(intial_cond_soliton(x,c,0));
B_0 = g.*fft(real(ifft(phi_hat_0)).^2);

t = h; step = 1;
%set(gcf,'doublebuffer','on')
%v = VideoWriter('yj627_SBDF2.avi');
%open(v);
while t + h/2 < T
    step = step+1;
    t = t+h;
    B_1 = g.*fft(real(ifft(phi_hat_1)).^2);
    phi_tmp = ((4/3)*phi_hat_1 - (1/3)*phi_hat_0 + (2*h/3)*(2*B_1 - B_0))./invA; % SBDF2 step
    phi_hat_0 = phi_hat_1;
    phi_hat_1 = phi_tmp;
    B_0 = B_1;
    
    %if mod(step,1000)==0 % Plot every 1000 steps
        %u = M * real(ifft([phi_hat_1(1:(N-1)/2)./N; zeros(M-N,1); phi_hat_1((N+1)/2:N)./N])); set(p,'ydata',u)
        %title(sprintf('t = %7.5f',t),'fontsize',18), drawnow
        %writeVideo(v,getframe);
    %end
    
    %if t + h/2 > T
        %u = M * real(ifft([phi_hat_1(1:(N-1)/2)./N; zeros(M-N,1); phi_hat_1((N+1)/2:N)./N]));
        %plot(x_fine,u,'.-');
    %end
end

%close(v);

philis = [philis, phi_hat_1];
%err1 = [err1, e1_rel_err(f_fine,u,L)];
%err2 = [err2, e2_rel_err(f_fine,u,L)];
%err3 = [err3, inf_rel_err(f_fine,u)];
end

%loglog(hlis,err1,'o-'); hold on
%loglog(hlis,err2,'o-');
%loglog(hlis,err3,'o-')

errs = [];
for i=2:length(hlis)
    err = norm(philis(:,i)-philis(:,i-1));
    errs = [errs,err];
end