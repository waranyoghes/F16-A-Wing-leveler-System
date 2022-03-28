W=[11 0 0;0 0.01 0;0 0 11];
V=sqrt(100);

x0 = [3;-2;1];
x_hat0=[-1;-2;0];%initial condition.
P0=[1 0.5 0.9; 0.3 0.7 0.8;1 0.6 0.2];
Ts = 0.5;%sample period
h = 0.01; %small time step to propagate states
Tfinal = 25; %total time
nframes = Tfinal/h; 

%defining continuous-time system in state space
A = [0, -18 162; 1,0,0;0,1,-10];
B = [1;0; 0];
C = [0,18,-162];
D = 0;
G=eye(3);

% weighting matrices
sysc = ss(A,B,C,D);


%calculating discrete-time system with zero-order hold (zoh is the default)
sysd = c2d(sysc,0.01,'zoh');


%xk = zeros(3,nframes+1);
x_hat= zeros(3,nframes+1);
y= zeros(1,nframes+1);
x = zeros(3,nframes+1);
Pk_0=P0;
xk_0=x0;
for i = 0:nframes 

    if i == 0     

        x(:,i+1) = sysd.A*x0+G*W*randn(3,1);
        y(:,i+1)= sysd.C*x(:,i+1)+V*randn(1,1);
        
        Pk_=sysd.A*Pk_0*transpose(sysd.A)+G*W*transpose(G);
        xk=sysd.A*xk_0;
        
        Pk=Pk_-Pk_*transpose(sysd.C)*inv([sysd.C*Pk_*transpose(sysd.C)+V])*sysd.C*Pk_;
        x_hat(:,i+1) = xk + Pk*transpose(sysd.C)*inv(V)*[y(:,i+1)-sysd.C*xk];
            

        
    else        
        x(:,i+1) = sysd.A*x(:,i)+G*W*randn(3,1);
        y(:,i+1)= sysd.C*x(:,i+1)+V*randn(1,1);
        
        Pk_=sysd.A*Pk_0*transpose(sysd.A)+G*W*transpose(G);
        xk=sysd.A*xk_0;
        
        Pk=Pk_-Pk_*transpose(sysd.C)*inv([sysd.C*Pk_*transpose(sysd.C)+V])*sysd.C*Pk_;
        x_hat(:,i+1) = xk + Pk*transpose(sysd.C)*inv(V)*[y(:,i+1)-sysd.C*xk];
        

    end
    Pk_0=Pk;
    xk_0=x_hat(:,i+1);
end

hold on

figure(1)

plot(h*[0:nframes],[x0(1), x(1,1:end-1)],':','LineWidth',2)
plot(h*[0:nframes],[x_hat0(1), x_hat(1,1:end-1)],'--','LineWidth',2)

xlabel('time (s)')
grid on
%xlim([0,50])
ylabel('states1')
hold off

figure(2)

hold on
%stairs(h*[0:nframes],[u0,u(1:end-1)],'LineWidth',2)
plot(h*[0:nframes],[x0(2), x(2,1:end-1)],':','LineWidth',2)
plot(h*[0:nframes],[x_hat0(2), x_hat(2,1:end-1)],'--','LineWidth',2)
xlabel('time (s)')
grid on
ylabel('states2')
hold off

figure(3)

hold on
%stairs(h*[0:nframes],[u0,u(1:end-1)],'LineWidth',2)
plot(h*[0:nframes],[x0(3), x(3,1:end-1)],':','LineWidth',2)
plot(h*[0:nframes],[x_hat0(3), x_hat(3,1:end-1)],'--','LineWidth',2)
xlabel('time (s)')
grid on
ylabel('states3')
hold off