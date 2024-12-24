%% Initializing the Variables
M = 1000; % Mass of the cart
m1 = 100; % Mass of the 1st pendulum
m2 = 100; % Mass of the 2nd pendulum
l1 = 20;  % Length of the 1st pendulum
l2 = 10;  % Length of the 2nd pendulum
g = 9.81; % Acceleration due to gravity

%% Defining the A, B and C matrices
A = [0 1 0 0 0 0;
    0 0 -m1*g/M 0 -m2*g/M 0;
    0 0 0 1 0 0;
    0 0 (-m1*g-M*g)/(M*l1) 0 -m2*g/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -m1*g/(M*l2) 0 (-m2*g-M*g)/(M*l2) 0];

B = [0;
    1/M;
    0;
    1/(M*l1);
    0;
    1/(M*l2)];

C1 = [1, 0, 0, 0, 0, 0]; % Observing only x(t)
C2 = [1, 0, 0, 0, 1, 0; 0, 0, 0, 0, 1, 0]; % Observing x(t) and theta_1(t)
C3 = [1, 0, 1, 0, 1, 0; 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 1, 0]; % Observing x(t), theta_1(t), theta_2(t)

D = 0;

%% Define noise covariances for Kalman filter
W = 0.01 * eye(size(A)); % Process noise
V1 = 0.01 * eye(size(C1, 1)); % Measurement noise for C1
V2 = 0.01 * eye(size(C2, 1)); % Measurement noise for C2
V3 = 0.01 * eye(size(C3, 1)); % Measurement noise for C3

%% Initial state
x_0 = [0; 0; 0.1; 0; 0.2; 0]; % Initial state for linear systems
x_0_nonlinear = [x_0; zeros(6, 1)]; % Initial state for nonlinear system

%% Design LQR Controller
Q = 1000 * eye(6); % State weights
R = 0.1; % Control input weight
K = lqr(A, B, Q, R);

%% Design Kalman Filter Gains
LC1 = lqe(A, W, C1, W, V1);
LC2 = lqe(A, W, C2, W, V2);
LC3 = lqe(A, W, C3, W, V3);

%% Combine LQR and Kalman Gains for C1
A_lqg_C1 = [A - B*K, zeros(size(A)); LC1*C1, A - LC1*C1];
B_lqg = [B; zeros(size(B))];
C_lqg_C1 = [C1, zeros(size(C1))];
sys_C1 = ss(A_lqg_C1, B_lqg, C_lqg_C1, D);

% Simulate Initial Response and Step Response for C1
figure;
initial(sys_C1, [x_0; zeros(size(x_0))]);
title('Initial Response for C1');
grid on;

figure;
step(sys_C1);
title('Step Response for C1');
grid on;

%% Combine LQR and Kalman Gains for C2
A_lqg_C2 = [A - B*K, zeros(size(A)); LC2*C2, A - LC2*C2];
C_lqg_C2 = [C2, zeros(size(C2))];
sys_C2 = ss(A_lqg_C2, B_lqg, C_lqg_C2, D);

% Simulate Initial Response and Step Response for C2
figure;
initial(sys_C2, [x_0; zeros(size(x_0))]);
title('Initial Response for C2');
grid on;

figure;
step(sys_C2);
title('Step Response for C2');
grid on;

%% Combine LQR and Kalman Gains for C3
A_lqg_C3 = [A - B*K, zeros(size(A)); LC3*C3, A - LC3*C3];
C_lqg_C3 = [C3, zeros(size(C3))];
sys_C3 = ss(A_lqg_C3, B_lqg, C_lqg_C3, D);

% Simulate Initial Response and Step Response for C3
figure;
initial(sys_C3, [x_0; zeros(size(x_0))]);
title('Initial Response for C3');
grid on;

figure;
step(sys_C3);
title('Step Response for C3');
grid on;

%% Define Nonlinear Dynamics Function
function dydt = nonLinear(t, y, A, C, LC, K, M, m1, m2, l1, l2, g)
    F = -K * y(1:6); % Control input
    dydt = zeros(12, 1);

    % Nonlinear system dynamics
    dydt(1) = y(2);
    dydt(2) = (F - (g/2) * (m1 * sin(2*y(3)) + m2 * sin(2*y(5))) ...
             - m1 * l1 * (y(4)^2) * sin(y(3)) - m2 * l2 * (y(6)^2) * sin(y(5))) ...
             / (M + m1 * (sin(y(3))^2) + m2 * (sin(y(5))^2));
    dydt(3) = y(4);
    dydt(4) = (dydt(2) * cos(y(3)) - g * sin(y(3))) / l1;
    dydt(5) = y(6);
    dydt(6) = (dydt(2) * cos(y(5)) - g * sin(y(5))) / l2;

    % Estimation error dynamics
    dydt(7:12) = (A - LC*C) * y(7:12);
end

%% Simulate Nonlinear System with C1
tspan = 0:0.1:1000;
[t, y] = ode45(@(t, y) nonLinear(t, y, A, C1, LC1, K, M, m1, m2, l1, l2, g), tspan, x_0_nonlinear);

% Plot Results
figure;
plot(t, y(:, 1:6));
title('Nonlinear System Response with C1');
xlabel('Time (s)');
ylabel('States');
grid on;
