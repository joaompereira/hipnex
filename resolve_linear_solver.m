function [linsolver, is_iterative] = resolve_linear_solver(linsolver)
    % Approximate Linear solver of (J(x) + s I) z = p
    %
    % The method has the following arguments:
    % - method: Use the method to solve the linear system. Can be either:
    %          - 'solve': for solving the system exactly using LU
    %          - 'lsqr': for using the LSQR method (iterative)

    if isstring(linsolver)||ischar(linsolver)
        % linsolver is a string describing the method to use
        switch linsolver
            case 'solve'
                is_iterative = false;
                linsolver = @solve_wrapper;
            case 'lsqr'
                is_iterative = true;
                linsolver = @lsqr_wrapper;
            otherwise
                error('NotImplementedError');
        end
    else
        % method is function provided by user
        % Assume it is iterative solver
        is_iterative = true;
    end
end


function x = solve_wrapper(J, s, p, ~, ~) 
    if isnumeric(J)
        x = (s * J + eye(length(p))) \ p;
    else
        warning('Using a function handle for J with exact solve is not recommended')
        I = eye(length(p));
        x = (s * J(I) + I) \ p;
    end

end

function [x, its] = lsqr_wrapper(J, s, p, tol, x0) 
    
    n = size(p, 1);

    if isnumeric(J)
        M = s * J + eye(n);
    else
        M = @(p, varargin) s * J(p, varargin{:}) + p;
    end
    
    if nargin < 5 || isempty(x0)
        [x, ~, its]  = lsqrSOL(n, n, M, p, 0, 0, 0, tol, 1e12, max(1e4, n), 0);
    else
        if isnumeric(J)
            Mx0 = M * x0; 
        else 
            Mx0 = M(x0, 1); 
        end
        [x, ~, its] = lsqrSOL(n, n, M, p - Mx0, 0, 0, 0, tol, 1e12, max(1e4, n), 0);
        x = x + x0;
    end


end