function [x, stats] = plain_npe(F, x, J, linsolver, varargin)
    % This function implements a plain version of Newton 
    % Proximal-Extragradient (NPE) algorithm for solving monotone 
    % variational inequalities.
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
    %    sigma_l: Optimization hyperparameter (see implementation of NPE)
    %    sigma_u: Optimization hyperparameter (see implementation of NPE)
    %        tol: Stopping tolerance for the algorithm: NPE stops if
    %             || F(x) || < tol. Defaults to 1e-6
    %    maxiter: Maximum number of iterations of NPE
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
    addOptional(p, 'sigma_l', []);
    addOptional(p, 'sigma_u', []);
    addOptional(p, 'tol', 1e-6);
    addOptional(p, 'maxiter', 10000);
    addOptional(p, 'verbosity', 2);
    parse(p, varargin{:});

    L = p.Results.L;
    hat_sigma = p.Results.hat_sigma;
    sigma_l = p.Results.sigma_l;
    sigma_u = p.Results.sigma_u;
    tol = p.Results.tol;
    maxiter = p.Results.maxiter;
    verbosity = p.Results.verbosity;

    %% Statistics log
    save_stats = (nargout > 1);

    if save_stats
        F_norms = zeros(1, 1 + maxiter);
        lambdas = zeros(1, 1 + maxiter);
        it_times = zeros(1, 1 + maxiter);
        fcn_timer = tic;
    end

    %% Define J using finite diferences if not provided
    if nargin < 3 || isempty(J)
        J = finite_difference_wrapper(F, x, p.Results.finite_diff_tol);
    end

    %% Define the linear solver
    % Here we define what will be the solver for the linear system
    %  (s J(x) + I) z = p, where s, x, and p are inputs
    % it can be an exact solver or iterative
    [linsolver, is_iterative] = resolve_linear_solver(linsolver);
    
    %% Define the hyperparameters of the algorithm
    assert(0 <= hat_sigma && hat_sigma < 1, 'hat_sigma is out of bounds')
    if isempty(sigma_u)
        sigma_u = .9 * (1 - hat_sigma);
    else
        assert(0 < sigma_u && hat_sigma+sigma_u < 1, 'sigma_u is out of bounds');
    end
    if isempty(sigma_l)
        sigma_l = .5 * sigma_u * (1 - hat_sigma) / (1 + hat_sigma);
    else
        assert(0 < sigma_l && sigma_l*(1+hat_sigma) <= sigma_u*(1-hat_sigma), ...
                'sigma_l is out of bounds');
    end

    alpha_m = 2*sigma_l/L;
    alpha_p = 2*sigma_u/L;
    
    %% First definitions (before iterating)
    Fx = F(x);
    normFx = norm(Fx);

    if normFx < tol
        return;
    end

    % Initial print
    if verbosity == 2
        fprintf('| iter | norm_F |\n');
    elseif verbosity == 3
        if is_iterative
            fprintf('| iter | norm_F | lambda | inner |\n');
        else
            fprintf('| iter | norm_F | lambda |\n');
        end
    end

    inner_iter = [];

    tl = [];
    tu = [];

    lambda = sqrt(2*sigma_l/(L*normFx));

    if save_stats
        lap = toc(fcn_timer);
        it_times(1) = lap;
        F_norms(1) = normFx;
        lambdas(1) = lambda;
    end
    
    it = 0;
    print_diagnostics
    total_inner = 0;
    F_evals = 1;
    DF_evals = 0;
    linsolve_count = 0;

    %% Start of the iteration
    for it = 1:maxiter
        
        %% Bracketing
        lambda = sqrt(2*sigma_l/(L*normFx));
        
        Jx = J(x); 
        DF_evals = DF_evals + 1;

        if is_iterative
            [dx, inner_iter] = linsolver(Jx, lambda, lambda * Fx, hat_sigma);
            total_inner = total_inner + inner_iter;
        else
            dx = linsolver(Jx, lambda, lambda * Fx, hat_sigma);
        end
        linsolve_count = linsolve_count + 1;
        
        print_diagnostics
        
        norm_dx = norm(dx);

        do_bissection = true;

        if lambda * norm_dx > alpha_p    
            tl = alpha_m / norm_dx;
            tu = lambda;
        elseif lambda * norm_dx < alpha_m
            tl = lambda;
            tu = alpha_p / norm_dx;
        else
            do_bissection = false;
        end

        %% Bissection
        while do_bissection
            lambda = sqrt(tu * tl);
            
            % Solve Linear system
            if is_iterative
                [dx, inner_iter] = linsolver(Jx, lambda, lambda * Fx, hat_sigma);
                total_inner = total_inner + inner_iter;
            else
                dx = linsolver(Jx, lambda, lambda * Fx, hat_sigma);
            end

            linsolve_count = linsolve_count + 1;
            norm_dx = norm(dx);

            print_diagnostics           

            if lambda * norm_dx > alpha_p    
                tu = lambda;
            elseif lambda * norm_dx < alpha_m
                tl = lambda;
            else
                do_bissection = false;
            end
        end
    
        %% NPE step
        tilde_x = x - dx;
        Ftx = F(tilde_x);
        F_evals = F_evals + 1;

        normFtx = norm(Ftx);
        if normFtx < tol
            normFx = normFtx;
            x = tilde_x;
            break
        end
        
        x = x - lambda * Ftx;
        Fx = F(x);
        F_evals = F_evals + 1;
        normFx = norm(Fx);

        if normFx < tol
            break;
        end

        if save_stats
            new_lap = toc(fcn_timer);
            it_times(it+1) = new_lap - lap;
            F_norms(it+1) = normFtx;
            lambdas(it+1) = lambda;
            lap = new_lap;
        end

        
    end

    if save_stats
        new_lap = toc(fcn_timer);
        it_times(it+1) = new_lap - lap;
        F_norms(it+1) = normFtx;
        lambdas(it+1) = lambda;
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

    if verbosity > 0
        fprintf('      n_iters: %d\n', it); 
        fprintf(' final_norm_F: %.2E\n', normFx);
        fprintf('linear_solves: %d\n', linsolve_count);
        fprintf('      F_evals: %d\n', F_evals);
        fprintf('     DF_evals: %d\n', DF_evals);
        if is_iterative
            fprintf('  total_inner: %d\n', total_inner);
        end
    end

    function print_diagnostics
        % Print diagnostics
        if verbosity == 2
            fprintf('| %5d | %.2E|\n', it, normFx);
        elseif verbosity == 3
            if is_iterative
                if it > 0
                    fprintf('|%5d |%.2E|%.2E| %5d |\n', it, normFx, lambda, inner_iter);
                else
                    fprintf('|%5d |%.2E|%.2E|       |\n', it, normFx, lambda);
                end
            else
            fprintf('|%5d |%.2E|%.2E|\n', it, normFx, lambda);
            end
        end


    end


            
end
          