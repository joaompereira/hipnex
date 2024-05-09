%% MATLAB script for testing the performance of HIPNEX and NPE
% This script tests the performance of HIPNEX and NPE on a quadratic
% minimax problem. The problem is defined as follows:
%
%   min_x max_y 0.5*x'*A*x + x'*B*y - 0.5*y'*C*y
%
% where A, B, and C are random matrices. The script compares the performance
% of HIPNEX and NPE using an iterative method for solving the intermediate problems, 
% versus solving these exactly, using a LU decomposition.

n1 = 2000; % Dimension of x
n2 = 2000; % Dimension of y

seed = 5;
rng(seed);

[hF, hJ, hsolver] = quadratic_min_max_setup(n1, n2);

%% Setting hyperparameters for HIPNEX and NPE
hyperparameters.L = 1e-4;              % Lipschitz constant of the Hessian
hyperparameters.tol = 1e-6;            % Tolerance for the stopping criterion
hyperparameters.maxiter = int32(1e4);  % Maximum number of iterations
hyperparameters.verbosity = 1;         % Verbosity level

% Setting the initial guess
x0 = randn(n1+n2, 1);

fprintf('\n:: Using iterative method ::\n');

% Setting hat_sigma for the iterative methods
hyperparameters.hat_sigma = .15;

fprintf(':: HIPNEX ::\n');
start = tic;
hipnex(hF, x0, hJ, hsolver, hyperparameters);
fprintf('Total Time: %.2f\n', toc(start));

fprintf(':: NPE ::\n');
start = tic;
plain_npe(hF, x0, hJ, hsolver, hyperparameters);
fprintf('Total Time: %.2f\n', toc(start));


fprintf('\n:: Using exact method (LU) ::\n');

% Setting hat_sigma for the exact methods
hyperparameters.hat_sigma = 1e-8;

fprintf(':: HIPNEX ::\n');
start = tic;
hipnex(hF, x0, hJ, 'solve', hyperparameters);
fprintf('Total Time: %.2f\n', toc(start));

fprintf(':: NPE ::\n');
start = tic;
plain_npe(hF, x0, hJ, 'solve', hyperparameters);
fprintf('Total Time: %.2f\n', toc(start));
