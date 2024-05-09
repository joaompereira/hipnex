function [hF, hJ, hJp, hsolver] = cubic_min_max_setup(n, L, kappa)
%% MATLAB function for setting up the cubic min-max problem
%
%   min_x max_y L * \|x\|^3 / 6 + y'*(A*x - b)
%
% where A is a square invertible matrix and b is a vector. Takes the following arguments:
%      n: dimension of A
%      L: Constant multiplying the cubic term, which coincides with the Lipschitz constant 
%         of the Hessian
%  kappa: Condition number of A
%
% It returns function handles, hF, hJ, hJp, and hsolver, which are the function, Jacobian, 
% Jacobian-vector product, and solver for the intermediate problem, respectively. The solver
% is only used when solving the inner problems iteratively, and uses the Jacobian-vector
% product, that consists of the function that calculates the multiplication by the Jacobian.

    if nargin < 2
        L = 0.1;
    end
    if nargin < 3
        kappa = 10;
    end

    log_kappa = log(kappa);

    [U, ~] = qr(randn(n, n));
    [V, ~] = qr(randn(n, n));
    S = diag(exp(linspace(-log_kappa, 0, n)));
    A = U * S * V;
    b = randn(n, 1) / sqrt(n);
    
    min_max_scale = [ones(n, 1); -ones(n, 1)];

    hF = @F;
    hJ = @J;
    hJp = @J_matvec;
    hsolver = @linsolver;

    % Defining function
    function F = F(z)
        x = z(1:n);
        y = z(n+1:end);
        nx = norm(x);
        F = [((L * nx / 2) * x + A' * y); (b - A * x)];
    end

    % Defining Jacobian
    % Reframing minmax as a variational inequality problem, implies that 
    % The Jacobian for the part that is being maximized is negated.
    function J = J(x)
        x1 = x(1:n);
        nx = norm(x1);

        upper_block = (L * nx / 2) * eye(n) + (L / (nx * 2)) * (x1 * x1');
        J = [upper_block, A'; -A, zeros(n, n)];
    end

    % Defining function that calculates J*p
    % where J is the Jacobian at x
    function J = J_matvec(x)
        x1 = x(1:n);
        nx = norm(x1);
        
        J = @matvec;

        function mv = matvec(z, t)
            z1 = z(1:n);
            z2 = z(n+1:end);
            if nargin < 2 || t==1
                mv_z1 = (L * nx / 2) * z1 + (L * dot(z1, x1) / (nx * 2)) * x1 + A' * z2;
                mv = [mv_z1; -(A * z1)];
            else
                mv_z1 = (L * nx / 2) * z1 + (L * dot(z1, x1) / (nx * 2)) * x1 - A' * z2;
                mv = [mv_z1; A * z1];
            end
        end
    end

    function [z, inner_iter] = linsolver(J, s, p, tol)
        Hxmz = @(z) min_max_scale .* (s*J(z) + z);
           
        mp = min_max_scale .* p;
        [z, ~, inner_iter] = minres_sol(Hxmz, mp, [], [], [], tol);
    end
    
end 