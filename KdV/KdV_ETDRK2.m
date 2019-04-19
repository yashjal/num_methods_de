N = 121; % # Fourier modes odd
M = N + 1000; % finer grid size
L = 60; % length of domain interval
c = 1; % speed of soliton
x = linspace(-L/2,L/2,N+1)'; x(end) = [];
x_fine = linspace(-L/2,L/2,M+1)'; x_fine(end) = [];
f_fine = intial_cond_soliton2(x_fine,c,3*c,0); %intial_cond_soliton(x_fine,c,0);
T = L/(c*2); % time period
h = 2e-4; % delta t

%hlis = 0.0001 * 2.^(0:10);
%err11 = []; err21 = []; err31 = [];
%philis = [];

ylimits=[-0.5 2];
figure(2); clf
p = plot(x_fine,f_fine,'linewidth',3);
xlim([-L/2 L/2])
ylim(ylimits);
grid on

% pre-assign values that doesn't change after each time step
k = (2*pi/L) * [0:(N-1)/2 -(N-1)/2:-1]';
A = 1i*k.^3;
g = -3i*k;
%for h=hlis
expA = exp(h*A);
predA = (expA-ones(N,1))./A; predA(1) = 0;
corrA = (expA-ones(N,1)-h*A)./(h*(A.^2)); corrA(1) = 0; 
phi_hat_0 = fft(intial_cond_soliton2(x,c,3*c,0));

t = 0; step = 0;
set(gcf,'doublebuffer','on')
v = VideoWriter('yj627_ETDRK2.avi');
open(v);
while t + h/2 < T
    step = step+1;
    t = t+h;
    B_0 = g.*fft(real(ifft(phi_hat_0)).^2);
    pred = expA.*phi_hat_0 + predA.*B_0; 
    B_pred = g.*fft(real(ifft(pred)).^2);
    corr = pred + corrA.*(B_pred-B_0); % ETDRK2 step
    phi_hat_0 = corr;
    
    if mod(step,1000)==0 % Plot every 1000 steps
        u = M * real(ifft([phi_hat_0(1:(N-1)/2)./N; zeros(M-N,1); phi_hat_0((N+1)/2:N)./N])); set(p,'ydata',u)
        title(sprintf('t = %7.5f',t),'fontsize',18), drawnow
        writeVideo(v,getframe);
    end
    
    %if t + h/2 > T
        %u = M * real(ifft([phi_hat_0(1:(N-1)/2)./N; zeros(M-N,1); phi_hat_0((N+1)/2:N)./N]));
        %plot(x_fine,u,'.-');
    %end    
end

close(v);

%philis = [philis, corr];
%err11 = [err11, e1_rel_err(f_fine,u,L)];
%err21 = [err21, e2_rel_err(f_fine,u,L)];
%err31 = [err31, inf_rel_err(f_fine,u)];
%end

%loglog(hlis,err11,'*-'); hold on
%loglog(hlis,err21,'*-');
%loglog(hlis,err31,'*-')

%errs1 = [];
%for i=2:length(hlis)
%    err = norm(philis(:,i)-philis(:,i-1));
%    errs1 = [errs1,err];
%end