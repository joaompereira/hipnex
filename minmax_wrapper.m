function [hF, hDF, hDFp, hlinsolver]= minmax_wrapper(n1, n2, grad, hess, hessp)
% Minmax solver wrapper for NPE and HIPNEX
% Converts a minmax problem to a variational inequality problem
% Also provides an iterative solver for the inner problems based on MINRES
% 
%  -- INPUTS --
%    n1: dimension of x (minimizing over x)
%    n2: dimension of y (maximizing over y)
%  grad: function that returns the gradient at (x,y). Should return  
%        a vector of length (n1+n2), with the blocks [g_x; g_y].
%  hess: function that returns the hessian at (x,y). Should return a
%        symmetric (n1+n2) x (n1+n2) matrix formed by the blocks
%                   | H_xx  H_xy |
%                   | H_xy  H_yy |
% hessp: function that returns the hessian vector product function at 
%        (x,y). That is, a function that takes a vector p = [p_x; p_y] as
%        input and returns H(x,y)p.
% 
% -- OUTPUTS --
%   hF: Handle for the variational inequality function.
%  hDF: Handle for the variational inequality Jacobian.
% hDFp: Handle for the variational inequality Jacobian vector product
%       function. Should take x as a vector and return DF*x
% hlinsolver: Linear solver for inner problems using MINRES

    min_max_scale = [ones(n1, 1); -ones(n2, 1)];
    hF = @(xy) min_max_scale .* grad(xy);
    if isempty(hess)
        hDF = [];
    else
        hDF = @(xy) min_max_scale .* hess(xy);
    end
    if nargin<4 || isempty(hessp)
        hDFp = hDF;
    else
        hDFp = @DFp_wrapper;
    end
    hlinsolver = @linsolver;

    function DFp = DFp_wrapper(xy)
        Hxy = hess(xy);
        DFp = @(p, varargin) min_max_scale .* Hxy(p, varargin{:});
    end

    function [z, inner_iter] = linsolver(J, s, p, tol)
        if isnumeric(J)
           Jxmz = @(z) min_max_scale .* (s*(J * z) + z);
        else
           Jxmz = @(z) min_max_scale .* (s*J(z) + z);
        end
           
        mp = min_max_scale .* p;
        [z, ~, inner_iter] = minres_sol(Jxmz, mp, [], [], [], tol);
    end

end