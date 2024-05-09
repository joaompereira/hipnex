%% MATLAB script for testing the performance of HIPNEX and NPE
% This script tests the performance of HIPNEX and NPE on a cubic
% minimax problem. The problem is defined as follows:
%
%   min_x max_y \|x\|^3 + y'*(A*x - b)
%
% where A is a square invertible matrix and b is a vector. The script compares the
% performance of HIPNEX and NPE using an iterative method for solving the intermediate
% problems, versus solving these exactly, using a LU decomposition.

clc
clearvars 

% increase verbosity for more info
verbosity = 1;

% Setting random seed for reproducibility
seed = 3;
rng(seed)

% Problem setup
L = 1e-3;                    % Lipschitz constant of the Hessian
kappa = 20;                  % Condition number of A
nvals = [1000, 2000, 5000];  % Dimensions (of A) to test 
nn = length(nvals);
tol = 1e-6;                  % Stopping Tolerance
n_iters = int32(1e4);


% Common hyperparameters
hyperparameters.L = L;               % Lipschitz constant of the Hessian
hyperparameters.tol = tol;           % Tolerance for the stopping criterion
hyperparameters.maxiter = n_iters;   % Maximum number of iterations
hyperparameters.verbosity = 0;       % Verbosity level


% Only hat_sigma changes for iterative vs exact
% hat_sigma control the error in solving the inner problems.
% This error is 0 when the inner problems are solved exactly, 
% for instance using an LU decomposition. Nevertheless, both
% HIPNEX and NPE allow for an error in solving the inner problems
hyperparameters_exact = hyperparameters;
hyperparameters_exact.hat_sigma = 1e-12;

hyperparameters_iterative = hyperparameters;
hyperparameters_iterative.hat_sigma = .15;

Algs = {'HIPNEX (iterative)', @(hF, x0, hJ, hJp, hsolver) ...
        hipnex(hF, x0, hJp, hsolver, hyperparameters_iterative);
        'NPE (iterative)', @(hF, x0, hJ, hJp, hsolver) ...
        plain_npe(hF, x0, hJp, hsolver, hyperparameters_iterative);
        'HIPNEX (exact)', @(hF, x0, hJ, hJp, hsolver) ...
        hipnex(hF, x0, hJ, 'solve', hyperparameters_exact);
        'NPE (exact)', @(hF, x0, hJ, hJp, hsolver) ...
        plain_npe(hF, x0, hJ, 'solve', hyperparameters_exact);
        'GO-2',  @(hF, x0, hJ, hJp, hsolver) ...
        ORN_ls_simple(hF, hJ, 1, 0, .5, .5, x0, n_iters, tol)
        };

nalgs = size(Algs, 1);

stats = cell(nalgs, nn);



fprintf('\n Numerical experiments for cubic min-max example\n')

for n_ind=1:nn

    n = nvals(n_ind);

    fprintf('\nn=%d\n', n)

    % Setup the cubic min-max problem
    [hF, hJ, hJp, hsolver] = cubic_min_max_setup(n, L, kappa);
    
    % Initial point initialization
    x0 = randn(2 * n, 1) / sqrt(n);
    
    for alg_ind = 1:nalgs

        fprintf(':: %s ::\n', Algs{alg_ind, 1});
        alg = Algs{alg_ind, 2};

        start = tic;
        [x, stats_it] = alg(hF, x0, hJ, hJp, hsolver);
        stats_it.total_time = toc(start);
        fprintf('Total Time: %f\n', stats_it.total_time);
        stats{alg_ind, n_ind} = stats_it;
    end
end

if ~exist('results', 'dir')
   mkdir('results')
end
save results/cubic_min_max_experiment.mat

