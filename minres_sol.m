function [x, r, it] = minres_sol(A_mv, b, x0, maxit, tol, tol_sol)
% minres - Minimum Residual Method for solving linear systems
%
% Syntax:
%   [x, r, it] = minres(A_mv, b, x0, maxit, tol, tol_sol)
%
% Inputs:
%   - A_mv: A matrix or a function handle that performs matrix-vector multiplication with A.
%   - b: Right-hand side vector of the linear system.
%   - x0: Initial guess for the solution (optional).
%   - maxit: Maximum number of iterations (optional).
%   - tol: Tolerance for the residual norm convergence criterion (optional).
%   - tol_sol: Tolerance for the solution norm convergence criterion (optional).
%
% Outputs:
%   - x: Solution vector.
%   - r: Residual vector.
%   - it: Number of iterations performed.
%
% Description:
%   The minres function implements the Minimum Residual Method for solving linear systems of the form Ax = b,
%   where A is a matrix or a function handle that performs matrix-vector multiplication with A.
%   The method iteratively improves the solution x by minimizing the residual norm ||b - Ax||.
%   The algorithm terminates when the residual norm is below the specified tolerance (tol) or the solution norm is below the specified tolerance (tol_sol).
%   If no initial guess (x0) is provided, the function initializes x with zeros.
%   If no maximum number of iterations (maxit) is provided, the function uses the length of b as the default value.
%   The function returns the solution vector x, the residual vector r, and the number of iterations performed (it).
%
% Example:
%   A = rand(100);
%   b = rand(100, 1);
%   [x, r, it] = minres(A, b);
%
% References:
%   - Paige, C. C., & Saunders, M. A. (1975). Solution of sparse indefinite systems of linear equations. SIAM Journal on Numerical Analysis, 12(4), 617-629.
%   - Freund, R. W. (1993). Conjugate gradient-type methods for linear systems with complex symmetric coefficient matrices. SIAM Journal on Scientific Computing, 14(6), 1372-1406.
    if nargin < 3
        x0 = [];
    end
    if nargin < 4
        maxit = [];
    end
    if nargin < 5
        tol = 0;
    end
    if nargin < 6
        tol_sol = 0;
    end

    if ~isa(A_mv, 'function_handle')
        A_mv = @(s) A_mv * s;
    end

    if isempty(x0)
        x = zeros(size(b));
        r = b;
    else
        x = x0;
        r = b - A_mv(x0);
    end

    p0 = r;
    s0 = A_mv(p0);
    ss1 = norm(s0);
    p1 = p0 / ss1;
    s1 = s0 / ss1;

    if isempty(maxit)
        maxit = length(b);
    end

    for it = 1:maxit
        ss1 = norm(s0);
        p2 = p1;
        p1 = p0 / ss1;
        s2 = s1;
        s1 = s0 / ss1;

        alpha = dot(r, s1);
        x = x + alpha * p1;
        r = r - alpha * s1;

        r_norm = norm(r);
        if r_norm < tol
            break;
        end
        if r_norm < tol_sol * norm(x)
            break;
        end

        p0 = s1;
        s0 = A_mv(s1);
        beta1 = dot(s0, s1);
        p0 = p0 - beta1 * p1;
        s0 = s0 - beta1 * s1;

        if it > 1
            beta2 = dot(s0, s2);
            p0 = p0 - beta2 * p2;
            s0 = s0 - beta2 * s2;
        end
    end
end