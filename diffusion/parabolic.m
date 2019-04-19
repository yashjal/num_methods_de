% 1.1.2 hw5 parabolic pde
% solving for u(x,t) and estimating order of accuracy

d = @(x) 2 + cos(pi*x); % d(x)
usol = @(x,t) exp(-pi^2*t/4)*cos(0.5*pi*x); % manufactured soln
fxt = @(x,t) exp(-pi^2*t/4) * (0.5*pi^2*cos(0.5*pi*x) + 0.5*pi*(-pi*sin(pi*x).*sin(0.5*pi*x) + 0.5*pi*cos(pi*x).*cos(0.5*pi*x))); % f(x,t)
ht = @(t) exp(-pi^2*t/4); % h(t)
gt = @(t) -0.5*pi*exp(-pi^2*t/4); % g(t)
%gt = @(t) 0; % dirichlet
x0 = 0; % starting point
T = 1/4;

dt = T/100;
N = 16;
dx = 1/N;

err_stan = [];
err_stag = [];

%for i=1:8
    %dt = dt/2;
    %N = 2*N;
    %%dx = 1/N;
    x_stan = dx:dx:1; % includes Neumann boundary value, N points
    %x_stan = 0:dx:1; % includes boundary values, N+1 points
    %x_stan = dx:dx:(1-dx); % dirichlet both ends
    x_stag = dx/2:dx:(1-dx/2); % doesn't include boundary values, N-1 points
    un_ex_st = usol(x_stan,0)';
    un_ex_sg = usol(x_stag,0)';
    un_im_st = un_ex_st;
    un_im_sg = un_ex_sg;
    fx = @(x)fxt(x,0);
    [Ast,fst] = discretize(0,ht(0),gt(0),N,dx,fx,d,x0); % standard grid
    [Asg,fsg] = discretize(1,ht(0),gt(0),N,dx,fx,d,x0); % staggered grid
    fst1 = fst;
    fsg1 = fsg;

    % set(gcf,'doublebuffer','on')
    % v = VideoWriter('stag_cn.avi');
    % open(v);

    figure(1)
    plot(x_stan,usol(x_stan,0))
    drawnow()
    %writeVideo(v,getframe)
    pause(0.1)

    for t = dt:dt:T

        %fst1(1) = -ht(t);
        %fst1(N+1) = fxt(1,t)+2*d(1+dx/2)*gt(t)/dx;
        %fst1(2:N) = fxt(x0+((2:N)'-1)*dx,t);

        fst1(2:N-1) = fxt(x0+(2:(N-1))'*dx,t);
        fst1(N) = fxt(1,t)+2*d(1+dx/2)*gt(t)/dx;
        %fst1(2:(N-2)) = fxt(x0+(2:(N-2))'*dx,t);
        %fst1(N-1) = fxt(1-dx,t)+d(1-dx/2)*gt(t)/dx^2;
        fst1(1) = fxt(x0+dx,t)+ht(t)*d(x0+dx/2)/(dx^2);

        fsg1(1) = fxt(x0+dx/2,t)+2*d(x0)*ht(t)/(dx^2);
        fsg1(N) = fxt(1-dx/2,t)+d(N*dx)*gt(t)/dx;
        fsg1(2:N-1) = fxt(x0+((2:N-1)'-1)*dx+dx/2,t);

        un_ex_st = fe(Ast,un_ex_st,fst,dt);
        un_ex_sg = fe(Asg,un_ex_sg,fsg,dt);

        un_im_st = cn(Ast,un_im_st,fst,fst1,dt);
        un_im_sg = cn(Asg,un_im_sg,fsg,fsg1,dt);
        %un_im_st = be(Ast,un_im_st,fst1,dt);
        %un_im_sg = be(Asg,un_im_sg,fsg1,dt);

        fst = fst1;
        fsg = fsg1;

        figure(1)
        plot(x_stan,usol(x_stan,t)); hold on
        plot(x_stan',un_im_st); hold off
        axis([0 1 0 1])
        drawnow()
        %writeVideo(v,getframe)
        pause(0.1)
        
    end
    %err_stan = [err_stan; norm(un_im_st-(usol(x_stan,T))')];
    %err_stag = [err_stag; norm(un_im_sg-(usol(x_stag,T))')];
%end
%close(v)
%figure(1)
%loglog(2.^(0:7)*16,err_stag,'ko-')
%figure(2)
%loglog(2.^(0:7)*16,err_stan,'ro-')
%plot(x_stan,un_ex_st,'o-'); hold on
%plot(x_stag,un_ex_sg,'o-')
%plot(x_stan,un_im_st,'.-'); hold on
%plot(x_stag,un_im_sg,'.-')
%plot(x_stan,usol(x_stan,T))
%plot(x_stag,usol(x_stag,T))
%xlabel('x')
%ylabel('u(x,T)')
%legend('FE stan','FE stag','CN stan','CN stag','actual soln stan','actual soln stag')