function [hF, hJ, hsolver] = quadratic_min_max_setup(n1, n2)
%% MATLAB function for setting up the quadratic min-max problem
%
%   min_x max_y 0.5*x'*A*x + x'*B*y - 0.5*y'*C*y
%
% where A, B, and C are random matrices. Takes two arguments, n1 and n2,
% which are the dimensions of the matrices A and C. It returns three
% function handles, hF, hJ, and hsolver, which are the function, Jacobian,
% and an iterative solver for the quadratic problem relying on MINRES.
    
    A = randn(n1, n1);
    A = A * A' / n1 + eye(n1);
    B = randn(n1, n2) / sqrt(max(n1, n2));
    C = randn(n2, n2);
    C = C * C' / -n2 - eye(n2);   

    H = [A, B; -B', -C];
    H_sym = [A, B; B', C];
    min_max_scale = [ones(n1, 1); -ones(n2, 1)];

    hF = @F;
    hJ = @J;
    hsolver = @solve_linear_approx;

    %% Function definitions
    
    function Fx = F(x)
        Fx = H*x;
    end

    function Jx = J(x)
        Jx = H;
    end

    function [z, inner_iter] = solve_linear_approx(H, s, p, tol)
        
        % matvec operator of 
        Hmz = @(z) s*(H_sym * z) + min_max_scale .* z;
        
        mp = min_max_scale .* p;
        [z, ~, inner_iter] = minres_sol(Hmz, mp, [], [], [], tol);
    end

    
end

