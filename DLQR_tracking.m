
x0 = [0;0]; %initial condition. 
Ts = 0.5;%sample period
h = 0.1; %small time step to propagate states
Tfinal = 30; %total time
nframes = Tfinal/h; 

%defining continuous-time system in state space
A = [0, 1; -1, -1];
B = [0; 1];
C = [1 0];
D = 0;

% weighting matrices
Q = [1 0; 0 1];
R = 1;
N=[1;1];
sysc = ss(A,B,C,D);


%calculating discrete-time system with zero-order hold (zoh is the default)
sysd = c2d(sysc,h,'zoh');

[K,S,e] = dlqr(A,B,Q,R); %SDR gain
ref=[1;0];
%[Q_hat, R_hat, M] = find_SDR_QR(sysc,Q,R,Ts);

[P,K_,L,info] = idare(A,B,Q,R,[],[]);
Z=-inv(eye(2)+transpose(A)*P*B*inv(R)*transpose(B))*transpose(A)*P*(ref+pinv(B*inv(R)*transpose(B))*ref);
K1=-inv(R)*transpose(B)*P*inv(eye(2)+B*inv(R)*transpose(B)*P)*A;
K2=inv(R)*transpose(B)*P*[pinv(eye(2)+B*inv(R)*transpose(B)*P)*(B*inv(R)*transpose(B)*Z+ref)-Z];
disp(size(Z))
disp(size(K1))
disp(size(K2))


sysd = c2d(sysc,h,'zoh'); %using short time step for propagating states during simulation
ref = [0;0]; %step (reference input)

u = zeros(1,nframes+1);
y = zeros(2,nframes);
x = zeros(2,nframes+1);

for i = 0:nframes 

    if i == 0     

        u0 = K1*x0+K2;

        x(:,i+1) = sysd.A*x0 + sysd.B*u0;

        y0 = sysd.C*x0 + sysd.D*u0;

    else        

        if mod(i*h,Ts) == 0 %zero order hold. 

            u(i) = K1*x(:,i)+ K2; %feedback equation. 

        else

            if i==1

                u(i) = u0;

            else

                u(i) = u(i-1);

            end
        end

        x(:,i+1) = sysd.A*x(:,i) + sysd.B*u(i); %state update. 
        y(:,i) = sysd.C*x(:,i) + sysd.D*u(i);

    end
end

hold on

figure(1)

plot(h*[0:nframes],[x0(1), x(1,1:end-1)],'LineWidth',2)
plot(h*[0:nframes],[x0(2), x(2,1:end-1)],'LineWidth',2)

xlabel('time (s)')
grid on
%ylim([0,1.6])
ylabel('phi')
hold off

figure(2)
%stairs(h*[0:nframes],[u0,u(1:end-1)],'LineWidth',2)
plot(h*[0:nframes],[u0,u(1:end-1)],'LineWidth',2)
%ylim([-0.4,1.2])
xlabel('time (s)')
grid on
ylabel('control')



