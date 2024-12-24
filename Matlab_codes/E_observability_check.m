% Observability check
% symbolic variables are defined here for solving the problem
syms m1 m2 l1 l2 g  M  real 

% A matrix from state space representation
A = [0, 1, 0, 0, 0, 0;
     0, 0, -m1*g/M, 0, -m2*g/M, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, (-m1*g-M*g)/(M*l1), 0, -m2*g/(M*l1), 0;
     0, 0, 0, 0, 0, 1;
     0, 0, -m1*g/(M*l2), 0, (-m2*g-M*g)/(M*l2), 0];

disp('Matrix A:');
disp(A);



% Number of states
n = size(A, 1); 
fprintf('Number of states: %d\n', n);

function O = observability_matrix(A, C, n)
% Function to find the observability matrix.
%
% The observability matrix is  O = [C; C*A; C*A^2; ...; C*A^(n-1)].
%
% Inputs:
%   A - (n x n matrix) State matrix from the state-space representation.
%   C - (p x n matrix) Output matrix from the state-space representation.
%   n - (integer) Number of states.
%
% Outputs:
%   O - (p*n x n matrix) Observability matrix.

    O = C; % Initializing with first row
    for i = 1:n-1
        O = [O; C * A^i]; 
    end
end



function check_observability(A, C, n)
%Function to check observability
    Observability_Matrix = observability_matrix(A, C, n);
    % Display the observability matrix, commented out to make display
    % cleaner uncomment to see the observability matrix.
    % disp('Observability Matrix:');
    % disp(Observability_Matrix);
    
    % Compute the rank of the observability matrix
    rank_O = rank(Observability_Matrix);
    fprintf('Rank of the Observability Matrix: %d\n', rank_O);
    
    % Check observability
    if rank_O == n
        disp('The system is observable.');
    else
        disp('The system is not observable.');
    end
end


disp('--------------------------------------------------');

% When output vector = x(t)
% C matrix from state space representation
C = [1,0,0,0,0,0];
disp('When output vector = x(t) :')
disp('Matrix C:');
disp(C);
check_observability(A, C, n);

disp('--------------------------------------------------');

% When output vector = (θ1(t), θ2(t))
% C matrix from state space representation
C = [0, 0, 1, 0, 0, 0;  % θ1(t)
     0, 0, 0, 0, 1, 0]; % θ2(t)

disp('When output vector = (θ1(t), θ2(t)) :')
disp('Matrix C:');
disp(C);
check_observability(A, C, n);

disp('--------------------------------------------------');


% When output vector = (x(t), θ2(t))
% C matrix from state space representation
C = [1, 0, 0, 0, 0, 0;  % x(t)
     0, 0, 0, 0, 1, 0]; % θ2(t)
disp('When output vector = (x(t), θ2(t)) :')
disp('Matrix C:');
disp(C);
check_observability(A, C, n);

disp('--------------------------------------------------');


% When output vector = (x(t), θ1(t), θ2(t))
% C matrix from state space representation
C = [1, 0, 0, 0, 0, 0;  % x(t)
     0, 0, 1, 0, 0, 0;  % θ1(t)
     0, 0, 0, 0, 1, 0]; % θ2(t)

disp('When output vector = (x(t), θ1(t), θ2(t)) :')
disp('Matrix C:');
disp(C);
check_observability(A, C, n);

disp('--------------------------------------------------');
