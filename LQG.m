W=10^-2*[1 0 0 0;0 0 1 0;0 0 1 0;0 0 0 0];
V=sqrt(10^-4)*[1 0; 0 1];

x0 = [1;-1;1;1];
x_hat0=[-2;3;2;4];%initial condition.
P0=10^-1*[1 0 0 0;0 0 1 0;0 0 1 0;0 0 0 0];
Ts = 0.1;%sample period
h = 0.05; %small time step to propagate states
Tfinal = 25; %total time
nframes = Tfinal/h; 

%defining continuous-time system in state space
A = [-0.128 0.0803 -0.997 0.035 ;-62.867 -1.306 0.297 0 ;8.277 0.00472 -0.259 0 ;0 1 0.0806 0 ];
B = [0.00746 0.0344 0 0.0143;-37.467 35.645 9.425 2.45;0.102 3.221 0.654 -2.437; 0 0 0 0];
C = [0 0 1 0;0 0 0 1];
D = 0;
G=eye(4);

% weighting matrices
sysc = ss(A,B,C,D);


%calculating discrete-time system with zero-order hold (zoh is the default)
sysd = c2d(sysc,Ts,'zoh');
%%%
Q = [300 0 0 0;0 5 0 0;0 0 0 0;0 0 0 100];
R = [100 0 0 0;0 400 0 0;0 0 100 0;0 0 0 150];

[K,S,e] = dlqr(sysd.A,sysd.B,Q,R);

%xk = zeros(3,nframes+1);
u = zeros(4,nframes+1);
ref=[0;0;0;0];
x_hat= zeros(4,nframes+1);
y= zeros(2,nframes+1);
x = zeros(4,nframes+1);
Pk_0=P0;
xk_0=x0;

for i = 0:nframes 

    if i == 0     
        u0=ref;
        %u0 = K*(ref - x_hat0);
        x(:,i+1) = sysd.A*x0+ sysd.B*u0 + G*W*randn(4,1);
        y(:,i+1)= sysd.C*x(:,i+1)+V*randn(2,1);
        
        Pk_=sysd.A*Pk_0*transpose(sysd.A)+G*W*transpose(G);
        xk=sysd.A*xk_0+sysd.B*u0;
        
        Pk=Pk_-Pk_*transpose(sysd.C)*inv([sysd.C*Pk_*transpose(sysd.C)+V])*sysd.C*Pk_;
        x_hat(:,i+1) = xk + Pk*transpose(sysd.C)*inv(V)*[y(:,i+1)-sysd.C*xk];
            

        
    else
        u(:,i)=ref;
        %u(:,i)=K*(ref-x_hat(:,i));
        x(:,i+1) = sysd.A*x(:,i)+ sysd.B*u(:,i) +G*W*randn(4,1);
        y(:,i+1)= sysd.C*x(:,i+1)+V*randn(2,1);
        
        Pk_=sysd.A*Pk_0*transpose(sysd.A)+G*W*transpose(G);
        xk=sysd.A*xk_0+sysd.B*u(:,i);
        
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

figure(4)

hold on
%stairs(h*[0:nframes],[u0,u(1:end-1)],'LineWidth',2)
plot(h*[0:nframes],[x0(4), x(4,1:end-1)],':','LineWidth',2)
plot(h*[0:nframes],[x_hat0(4), x_hat(4,1:end-1)],'--','LineWidth',2)
xlabel('time (s)')
grid on
ylabel('states3')
hold off