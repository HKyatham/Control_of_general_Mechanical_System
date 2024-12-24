%Luenberger design and simulations for output vector (x(t), θ1(t) ,θ2(t))
% Parameters from section D
M = 1000;    % Mass of the cart
m1 = 100;    % Mass of pendulum 1
m2 = 100;    % Mass of pendulum 2
l1 = 20;     % Length of pendulum 1
l2 = 10;     % Length of pendulum 2
g = 9.81;    % Gravitational acceleration

% A matrix from state space representation
A = [0, 1, 0, 0, 0, 0;
     0, 0, -m1*g/M, 0, -m2*g/M, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, (-m1*g-M*g)/(M*l1), 0, -m2*g/(M*l1), 0;
     0, 0, 0, 0, 0, 1;
     0, 0, -m1*g/(M*l2), 0, (-m2*g-M*g)/(M*l2), 0];

% B matrix from state space representation
B = [0;
     1/M;
     0;
     1/(M*l1);
     0;
     1/(M*l2)];



% Output: (x(t), θ1(t) ,θ2(t))
C = [1, 0, 0, 0, 0, 0;    % x(t)
      0, 0, 1, 0, 0, 0;    % θ1(t)
      0, 0, 0, 0, 1, 0];   % θ2(t)


D = 0;  

% Settings for simulation
t = 0:0.01:10;         % Time vector deom 0 to 10 seconds
u = ones(size(t));     % Unit step input
x0 = [0.1; 0 ; 0.1; 0; 0.1; 0];  % Initial condition for states
x_hat0 = [0; 0; 0; 0; 0; 0];    % Initial condition for observer
x0_nonlinear = [x0; zeros(6, 1)]; % Initial state for nonlinear system

tolerance = 1e-6;

% Objective function for optimization, it finds the settling time and
% sum of estimation error.
function [settling_time, sum_estimation_error] = objective(poles, A, B, C, u, t, x0, x_hat0, tolerance)
    % Extracting poles for C 
    poles_observer_C = poles;

    % Computing observer gains using pole placement method in mATLAB
    L1 = place(A', C', poles_observer_C)';
    
    % Simulating the actual system
    system_actual = ss(A, B, eye(6), zeros(6, 1));
    [x_actual, ~, x_states] = lsim(system_actual, u, t, x0);
    
    % Augmenting the matrices
    A_augmented = A - L1*C;
    B_augmented = [B, L1];
    % Defining the combined input matrix for observer system (including x, θ1, θ2)
    input_observer = [u(:), x_actual(:, 1), x_actual(:, 3), x_actual(:, 5)];  
    system_observer = ss(A_augmented, B_augmented, eye(6), zeros(6, 4));
    [~, ~, x_hat_states] = lsim(system_observer, input_observer, t, x_hat0); 
   
   
    % Calculating the estimation error
    error = x_states - x_hat_states;


    % Settling time calculation: find the when error is within tolerance of final value
    final_error = error(end, :);
    settling_times = zeros(1, 6);
    for i = 1:6
        error_each_state = error(:, i);
        threshold = abs(final_error(i)) * (1 + tolerance);  % Error tolerance
        
        % Find when the error first comes and stays within tolerance
        index_settle = find(abs(error_each_state) < threshold, 1, 'first');
        
        if ~isempty(index_settle)
            % Finding when the error stays within tolerance
            index_end = find(abs(error_each_state(index_settle:end)) > threshold, 1, 'first') + index_settle - 1;
            if isempty(index_end)  % If error stays within tolerance till end
                index_end = length(error_each_state);
            end
            settling_times(i) = t(index_end) - t(index_settle);  % Settling time in seconds
        end
    end
    settling_time = mean(settling_times);  % Average settling time across states

    % Sum of the estimation error values until settling
    sum_estimation_error = sum(sum(abs(error(1:index_end, :))));  
end

% Sum error function for optimization (returns sum of estimation error)
function sum_estimation_error = sum_error(poles, A, B, C, u, t, x0, x_hat0, tolerance)
    [~, sum_estimation_error] = objective(poles, A, B, C, u, t, x0, x_hat0, tolerance);
end

% Randomly chosen initial poles
initial_poles = [-1, -2, -3, -4, -5, -6];

w1 = 1;  % Weight for settling time
w2 = 0.5;  % Weight for sum of error



function penalty = pole_penalty(poles)
    penalty = 0;
    for i = 1:length(poles)
        if poles(i) < -1000
            penalty = penalty + (poles(i) + 1000)^2; % Penalty for being less than -1000
        elseif poles(i) > -0.0000
            penalty = penalty + (poles(i) + 0.1)^2; % Penalty for being greater than 0
        end
    end
end

loss_function = @(poles) w1 * objective(poles, A, B, C, u, t, x0, x_hat0, tolerance) + w2 * sum_error(poles, A, B, C, u, t, x0, x_hat0, tolerance) + pole_penalty(poles);  

% Optimization settings 
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 10000, 'MaxIter', 5000);
optimized_poles = fminsearch(loss_function, initial_poles, options);



% Displaying the optimized poles
disp('Optimized poles:');
disp(optimized_poles);
% Simulating and plottting the observer's performance with the optimized poles
poles_observer_C = optimized_poles(1:6);
L1 = place(A', C', poles_observer_C)';
% Displaying the "best" Lueneberger matrix
disp('Best observer gain matrix , L:');
disp(L1)




%% Simulating observer 
A_augmented = A - L1*C;
B_augmented = [B, L1];  % Augmenting B matrix to include observer feedback

% Simulating the actual system
system_actual_C = ss(A, B, eye(6), zeros(6, 1));
[x_actual_C, ~, x_states_C] = lsim(system_actual_C, u, t, x0);

% Defining the combined input matrix for observer system (including x, θ1, θ2)
input_observer = [u(:), x_actual_C(:, 1), x_actual_C(:, 3), x_actual_C(:, 5)];  % Augmented with x, θ1, θ2
% Simulating observer system for C2
sys_observer = ss(A_augmented, B_augmented, eye(6), zeros(6, 4));
[~, ~, x_hat_states] = lsim(sys_observer, input_observer, t, x_hat0);


% Calculating the estimation error
error = x_states_C - x_hat_states;
Q = 1000*eye(6);
R = 0.1;
K = lqr(system_actual_C,Q,R);

% Plotting 
figure;
set(gcf, 'Position', [100, 100, 1000, 1600]);  

for i = 1:6
    subplot(6, 1, i);  
    plot(t, error(:, i), 'LineWidth', 1.5);
    xlabel('Time (s)');
    if i == 1
        ylabel('Error in x(t)', 'Interpreter', 'latex');
    elseif i == 2
        ylabel('Error in $\dot{x}(t)$', 'Interpreter', 'latex');  
    elseif i == 3
        ylabel('Error in $\theta_1(t)$', 'Interpreter', 'latex');
    elseif i == 4
        ylabel('Error in $\dot{\theta}_1(t)$', 'Interpreter', 'latex'); 
    elseif i == 5
        ylabel('Error in $\theta_2(t)$', 'Interpreter', 'latex');
    elseif i == 6
        ylabel('Error in $\dot{\theta}_2(t)$', 'Interpreter', 'latex'); 
    end
    grid on;
end
sgtitle('Estimation Error for output vector $(x(t),\theta_1(t),\theta_2(t))$', 'Interpreter', 'latex');



%% Define Nonlinear Dynamics Function
function dydt = nonLinear(~, y, K, A, C, LC, M, m1, m2, l1, l2, g)
    F = -K*y(1:6); % Control input
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
[t, y] = ode45(@(t, y) nonLinear(t, y, K, A, C, L1, M, m1, m2, l1, l2, g), tspan, x0_nonlinear);

% Plot Results
figure;
plot(t, y(:, 1:6));
sgtitle('Nonlinear System Response with output vector $(x(t),\theta_1(t),\theta_2(t))$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('States');
grid on;