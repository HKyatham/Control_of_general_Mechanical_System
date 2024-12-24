%% Initializing the Variables
M = 1000; % Mass of the cart
m1 = 100; % Mass of the 1st pendulum
m2 = 100; % Mass of the 2nd pendulum
l1 = 20;  % Length of the 1st pendulum
l2 = 10;  % Length of the 2nd pendulum
g = 9.81; % Acceleration due to gravity

%% Defining the A, B and C matrices.
A = [0 1 0 0 0 0;
    0 0 -m1*g/M 0 -m2*g/M 0;
    0 0 0 1 0 0;
    0 0 (-m1*g-M*g)/(M*l1) 0 -m2*g/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -m1*g/(M*l2) 0 (-m2*g - M*g)/(M*l2) 0];

B = [0;
    1/M;
    0;
    1/(M*l1);
    0;
    1/(M*l2)];
C = eye(6);
D = 0;
% Controllability Matrix
CtrlM = ctrb(A,B);
% Rank of the controllability matrix
if (rank(CtrlM) == 6)
    disp("Rank of the controllability matrix is equal to n, System " + ...
        "is Controllable")
else
    disp("Rank of the controllability matrix is less than n, System " + ...
        "is not Controllable")
end
%% Defining Q and R values
Q = 1000*eye(6);
R = 0.1;
x_0 = [0; 0; 0.1; 0; 0.2; 0]; % Initial condition as a column vector
% Defining the system
system = ss(A,B,C,D);
% Calculating K, Solution to Riccati equation and Poles
[K,Solution,Poles] = lqr(system,Q,R);



figure;
initial(system, x_0); % Simulate states
grid on;

% Closed-loop system dynamics
Acl = A - B*K; % Closed-loop A matrix
system_cl = ss(Acl, B, C, D);

% Use initial() to simulate states
tspan = 0:0.01:30000; % Time vector
[y, t, x] = initial(system_cl, x_0, tspan); % Simulate states

% Compute the control input u(t)
u = -K * x'; % K * x' gives u for all time steps

% Plot states and control input
figure;
initial(system_cl, x_0)
grid on;
% Plot each state as a subplot
% for i = 1:6
%     subplot(7, 1, i); % 7 rows, 1 column, current subplot index
%     plot(t, x(:, i));
%     title(['State x_', num2str(i)]);
%     xlabel('Time (s)');
%     ylabel(['x_', num2str(i)]);
%     grid on;
% end
% 
% % Plot control input
% subplot(7, 1, 7);
% plot(t, u);
% title('Control Input u(t)');
% xlabel('Time (s)');
% ylabel('u(t)');
% grid on;
%% Define Nonlinear System Dynamics with LQR Control
function dydt = nonLinear(t, y, K, M, m1, m2, l1, l2, g)
    % LQR Feedback Control
    F = -K * y; % Control input based on current state
    
    % Nonlinear dynamics of the system
    dydt = zeros(6,1);
    dydt(1) = y(2);
    %y(2)=xdot;
    dydt(2)=(F-(g/2)*(m1*sind(2*y(3))+m2*sind(2*y(5)))...
        -(m1*l1*(y(4)^2)*sind(y(3)))-(m2*l2*(y(6)^2)*sind(y(5))))...
        /(M+m1*((sind(y(3)))^2)+m2*((sind(y(5)))^2)); %X_DD
    %y(3)=theta1;
    dydt(3)= y(4);
    %y(4)=theta1dot;
    dydt(4)= (dydt(2)*cosd(y(3))-g*(sind(y(3))))/l1'; %theta 1 Ddot;
    %y(5)=theta2;
    dydt(5)= y(6); %theta 2D
    %y(6)=theta2dot;
    dydt(6)= (dydt(2)*cosd(y(5))-g*(sind(y(5))))/l2; %theta 2Ddot;
end
figure;
%% Simulate the Nonlinear System Using ODE45
tspan = 0:0.01:10000; % Time span for simulation
[t, y] = ode45(@(t, y) nonLinear(t, y, K, M, m1, m2, l1, l2, g), ...
    tspan, x_0);
%plotting the function output on a 2D graph
plot(t,y)
grid on
disp(eig(Acl))