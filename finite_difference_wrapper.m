function J = finite_difference_wrapper(F, x, h)
    if nargin<3
        h = 1e-6;
    end
    Fx0 = F(x);
    J = @finite_difference;
    
    function [df] = finite_difference(p)
        
        df = (F(x0 + h*p) - Fx0) / h;
    end
end



