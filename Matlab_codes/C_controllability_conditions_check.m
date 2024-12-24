% Controllability conditions check
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

% B matrix from state space representation
B = [0;
     1/M;
     0;
     1/(M*l1);
     0;
     1/(M*l2)];

disp('Matrix B:');
disp(B);


% Number of states
n = size(A, 1); 
disp('Number of states:');
disp(n);

% 
function C = controllability_matrix(A, B, n)
% Function to find controllability Matrix
% The controllability matrix  = [B, AB, A^2B, ..., A^(n-1)B].
%
% Inputs:
%   A - (n x n matrix) A Matrix from state space representation
%   B - (n x m matrix) B Matrix from state space representation
%   n - (integer) Number of states.
%
% Outputs:
%   C - (n x n*m matrix) Controllability matrix.
    C = B; % Initialize
    for i = 1:n-1
        C = [C, A^i * B];
    end
end


Original_Controllability_Matrix = controllability_matrix(A, B, n);
% simplifying the expression. The simplifies expression is given in the
% report.
Simplified_Controllability_Matrix = simplify(Original_Controllability_Matrix);
disp('Controllability Matrix:');
disp(Simplified_Controllability_Matrix);

% Computing the rank of controllability matrix symbolically
rank_C = rank(Simplified_Controllability_Matrix);
fprintf('Symbolic Rank of the Controllability Matrix: %d\n', rank_C);

% Lets check the rank for specific conditions

% When the lengths of the pendulums are equal
% Applying substitution for l1 == l2
A_same_length = A;
B_same_length = B;
A_same_length = subs(A_same_length, l1, l2);  % Replace l1 with l2
B_same_length = subs(B_same_length, l1, l2);  % Replace l1 with l2

% When the mass of the pendulums are equal
% Applying substitution for m1 == m2
A_same_pendulum_mass = A;
B_same_pendulum_mass = B;
A_same_pendulum_mass = subs(A_same_pendulum_mass, m1, m2);  % Replace m1 with m2
B_same_pendulum_mass = subs(B_same_pendulum_mass, m1, m2);  % Replace m1 with m2


% When the mass of cart and mass of pendulums are same
% Applying substitution for m1 == m2 == M
A_same_all_mass = A_same_pendulum_mass;
B_same_all_mass = B_same_pendulum_mass;
A_same_all_mass = subs(A_same_all_mass, M, m2);  % Replace M with m2
B_same_all_mass = subs(B_same_all_mass, M, m2);  % Replace M with m2

% Display the updated matrices   % Used for coding ad debugging. Commenting
% here
% disp('Updated A_same_pendulum_mass matrix:');
% disp(A_same_pendulum_mass);
% disp('Updated B_same_pendulum_mass matrix:');
% disp(B_same_pendulum_mass);

% Display the updated matrices   
% disp('Updated A_same_all_mass matrix:');
% disp(A_same_all_mass);
% disp('Updated B_same_all_mass matrix:');
% disp(B_same_all_mass);

% Display the updated matrices
% disp('Updated A_same_length matrix:');
% disp(A_same_length);
% disp('Updated B_same_length matrix:');
% disp(B_same_length);

Controllability_Matrix = controllability_matrix(A_same_pendulum_mass, B_same_pendulum_mass, n);
% Computing the rank of controllability matrix symbolically
rank_C = rank(Controllability_Matrix);
fprintf('Rank of the Controllability Matrix when m1 = m2: %d\n', rank_C);


Controllability_Matrix = controllability_matrix(A_same_length, B_same_length, n);
% Computing the rank of controllability matrix symbolically
rank_C = rank(Controllability_Matrix);
fprintf('Rank of the Controllability Matrix when  l1 = l2: %d\n', rank_C);


Controllability_Matrix = controllability_matrix(A_same_all_mass, B_same_all_mass, n);
% Computing the rank of controllability matrix symbolically
rank_C = rank(Controllability_Matrix);
fprintf('Rank of the Controllability Matrix when  M = m1 = m2: %d\n', rank_C);
