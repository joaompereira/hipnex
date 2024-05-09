function [y, stats] = hipnex(F, x0, J, linsolver, varargin)
    % Implementation of 
    % A search-free O(1/k^(3/2)) Homotopy Inexact Proximal-Newton
    % EXtragradient algorithm for monotone variational inequalities
    % 
    % We focus on the unconstrained case, where, for a monotone operator F, 
    % we seek a solution to F(x) = 0;
    % 
    %% INPUT 
    % As a Newton-type method, it requires knowledge of the Jacobian of F.
    % The method has the following arguments:
    %         x0: Initial vector for optimization
    %          F: The monotone operator. Should be a function that given x
    %             as input, returns F(x)
    %          J: Function that calculates the Jacobian at x. Given x as 
    %             an input, it can return either a square matrix J, 
    %             (the Jacobian at x), or alternatively it can return 
    %             another function that accepts a vector p as an argument
    %             and calculates J p (this last version is particularly
    %             suited for iterative algorithms). For some iterative
    %             algorithms (such as LSQR, also included here) it is 
    %             convenient that this version of J also allows to
    %             calculate multiplication by the transpose (J^T p). This 
    %             should be done by adding an additional flag argument to
    %             the function, which should be 1 when calculating J p, 
    %             and 2 when calculating J^T p. If J is not provided, 
    %             finite-differences are used to approximate J p.
    %  linsolver: Function that calculates a solution to the shifted 
    %             equation (J + s I) z = p, where J is the Jacobian, s is 
    %             a scalar and p is a vector. It should take as inputs 
    %             J, s, p and tol, where J is the Jacobian, s is a scalar,
    %             p is a vector and tol is the tolerance for solving the 
    %             inner problem (as in eq (26) of the paper). It can also 
    %             accept strings as input, in that case it implements a 
    %             solver implemented in resolve_linear_solver. As default, 
    %             the inner problems are solved exactly
    % It also has the following optional arguments and hyper-parameters
    %          L: Lipschitz constant of the Jacobian. Defaults to 1.
    %  hat_sigma: Inner problem tolerance hyperparameter. Defauls to 0, as
    %             default solver is exact linear solver.
    %      theta: Algorithm hyperparameter introduced in equation (24)
    %      sigma: Algorithm hyperparameter introduced later in the
    %             analysis. The hyperparameter sigma is related to eta
    %             introduced in equation (25) by 
    %                   eta = 2 hat_theta / (sigma * L)
    %        tol: Stopping tolerance for the algorithm: HIPNEX stops if
    %             || F(x) || < tol. Defaults to 1e-6
    %    maxiter: Maximum number of iterations of HIPNEX
    %  verbosity: Verbosity of the algorithm. Is one of the following:
    %             0: Does not print any information
    %             1: Prints only information in the end of optimization
    %             2: Prints some information per iteration
    %             3: Prints more detailed information per iteration
    %             Defaults to 2.
    %
    %% OUTPUT 
    % The algorithm outputs the final y, where \|F(y)\| < tol, and some
    % statistics of the algorithm, saved in the stats struct. The
    % statistics are as follows:
    %    n_iters: Total number of iterations until convergence.
    % total_time: Total time until convergence.
    %   it_times: Computation time per iteration.
    %    F_norms: Norm of f at each iteration.
    %    lambdas: Value of lambda at each iteration.
    %    F_evals: Total number of function evaluations (F(x)).
    %   DF_evals: Total number of Jacobian evaluations (J(x)).
    % linsolve_count: Total number of linear problems solved.
    % total_inner: (Iterative methods only) Total number of inner
    %            iterations of iterative linear solvers.
    
    %% Setting up function parameters 
    p = inputParser;
    addOptional(p, 'L', 1);
    addOptional(p, 'hat_sigma', 0);
    addOptional(p, 'theta', []);
    addOptional(p, 'sigma', []);
    addOptional(p, 'tol', 1e-6);
    addOptional(p, 'maxiter', 10000);
    addOptional(p, 'verbosity', 2);
    parse(p, varargin{:});

    L = p.Results.L;
    hat_sigma = p.Results.hat_sigma;
    theta = p.Results.theta;
    sigma = p.Results.sigma;
    tol = p.Results.tol;
    maxiter = p.Results.maxiter;
    verbosity = p.Results.verbosity;
    
    %% Statistics log
    save_stats = (nargout > 1);

    if save_stats
        F_norms = zeros(1, maxiter+1);
        lambdas = zeros(1, maxiter+1);
        it_times = zeros(1, maxiter+1);
        fcn_timer = tic;
    end
    
    %% Define J using finite diferences if not provided
    if nargin < 3 || isempty(J)
        J = @(x) finite_difference_wrapper(F, x, p.Results.finite_diff_tol);
    end

    %% Define the linear solver
    % Here we define what will be the solver for the linear system
    %  (s J(x) + I) z = p, where s, x, and p are inputs
    % it can be an exact solver or iterative
    [linsolver, is_iterative] = resolve_linear_solver(linsolver);

    %% Define the hyperparameters of the algorithm
    % Assert (21) - (23) hold
    assert(0 <= hat_sigma && hat_sigma < .5, 'hat sigma violates (21)');
    if isempty(theta)
        theta = (1-hat_sigma)*(1-2*hat_sigma) / 2;
    else
        assert(0 < theta && theta < (1-hat_sigma)*(1-2*hat_sigma), 'theta violates (22)');
    end
    hat_theta = theta * (theta + hat_sigma - hat_sigma ^ 2) / (1 - hat_sigma) ^ 2;
    if isempty(sigma)
        sigma = .95;
    else
        assert(0 < sigma && sigma < 1, 'etaL violates (23)');
    end
    etaL = 2 * hat_theta / sigma;
    eta = etaL / L;
    tau = (theta - hat_theta) / (theta + etaL / 4 + sqrt((theta + etaL / 4) ^ 2 - theta*(theta-hat_theta)));
    
    %% First definitions (before iterating)
    x = x0;
    y = x;
    Fy = F(y);
    normFy = norm(Fy);

    if normFy < tol
        return;
    end

    % l is lambda in the paper
    lambda = sqrt(2*theta/(L*normFy));
    
    % Initial print
    if verbosity == 2
        fprintf('| iter | norm_F |\n');
        fprintf('|    0 |%.2E|\n', normFy);
    elseif verbosity == 3
        if is_iterative
            fprintf('| iter | norm_F | lambda | inner | HPE | lup |\n');
            fprintf('|    0 |%.2E|%.2E|       |     |     |\n', normFy, lambda);
        else
            fprintf('| iter | norm_F | lambda | HPE | lup |\n');
            fprintf('|    0 |%.2E|%.2E|     |     |\n', normFy, lambda);
        end
    end

    inner_iter = [];
    
    if save_stats
        lap = toc(fcn_timer);
        it_times(1) = lap;
        F_norms(1) = normFy;
        lambdas(1) = lambda;
    end

    total_hpe = 0;
    total_lup = 0;
    linsolve_count = 0;
    F_evals = 1;
    DF_evals = 0;
    total_inner = 0;
    
    %% Start of the iteration
    for it = 1:maxiter

        nhpe = 0;
        n_lambda_up = 0;

        % HIPNEX lambda's checks
        while norm(lambda * Fy + y - x) <= 2 * hat_theta / (lambda * L)

            while lambda * norm(y - x) >= eta
                x = x - (tau*lambda) * Fy;
                lambda = lambda * (1 - tau);
                nhpe = nhpe + 1;
            end

            lambda = lambda / (1 - tau);
            n_lambda_up = n_lambda_up + 1;

            % while l*L*norm(l*Fy+y-x) <= 2*theta
            %    l = l / (1-tau);
            %    n_lambda_up = n_lambda_up + 1;
            % end

        end

        total_lup = total_lup + n_lambda_up;
        total_hpe = total_hpe + nhpe;
        
        %% Newton iteration
        % Solve the regularized problem 
        r = lambda * Fy + y - x;
        Jy = J(y);
        DF_evals = DF_evals + 1; 

        if is_iterative
            [dy, inner_iter] = linsolver(Jy, lambda, r, hat_sigma);
            total_inner = total_inner + inner_iter;
        else
            dy = linsolver(Jy, lambda, r, hat_sigma);
        end
        linsolve_count = linsolve_count + 1;

        y = y - dy;
        Fy = F(y);
        F_evals = F_evals + 1; 

        normFy = norm(Fy);

        % Print diagnostics
        if verbosity == 2
            fprintf('| %5d | %.2E|\n', it, normFy);
        elseif verbosity == 3
            if is_iterative
                fprintf('|%5d |%.2E|%.2E| %5d |%4d |%4d |\n', it, normFy, lambda, inner_iter, nhpe, n_lambda_up);
            else
                fprintf('|%5d |%.2E|%.2E|%4d |%4d |\n', it, normFy, lambda, nhpe, n_lambda_up);
            end
        end
        
        % Save stats
        if save_stats
            new_lap = toc(fcn_timer);
            it_times(it+1) = new_lap - lap;
            F_norms(it+1) = normFy;
            lambdas(it+1) = lambda;
            lap = new_lap;
        end

        if normFy < tol
            break;
        end

    end

    if verbosity > 0
        fprintf('      n_iters: %d\n', it); 
        fprintf(' final_norm_F: %.2E\n', normFy);
        fprintf('linear_solves: %d\n', linsolve_count);
        fprintf('      F_evals: %d\n', F_evals);
        fprintf('     DF_evals: %d\n', DF_evals);
        if is_iterative
            fprintf('  total_inner: %d\n', total_inner);
        end
    end

    if save_stats
        stats.it_times = it_times(1:it+1);
        stats.F_norms = F_norms(1:it+1);
        stats.lambdas = lambdas(1:it+1);
        stats.linsolve_count = linsolve_count;
        stats.n_iters = it;
        stats.F_evals = F_evals;
        stats.DF_evals = DF_evals;
        if is_iterative
            stats.total_inner = total_inner;
        end
    end
end
