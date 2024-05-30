function varargout = mccormick(X)
% McCormick function 



    % if no input is given, return dimensions, bounds and minimum
    if (nargin == 0)
        varargout{1} = 2;  % # dims
        varargout{2} = [-1.5, -3]; % LB
        varargout{3} = [+4, +4]; % UB
        varargout{4} = [-5.471975602214493e-001, -1.547197559268372e+000]; % solution
        varargout{5} = -1.913222954981037e+000; % function value at solution

    % otherwise, output function value
    else
    
    % split input vector X into x1, x2        
    if size(X, 1) == 2
        x1 = X(1, :);        x2 = X(2, :); 
    else        
        x1 = X(:, 1);        x2 = X(:, 2);        
    end
    
    %   keep all values in the search domain
    x1(x1 < -1.5) = inf;       x1(x1 > 4) = inf;
    x1(x2 < -3.0) = inf;       x2(x2 > 4) = inf;
    
    % output function value
    varargout{1} = sin(x1 + x2) + (x1-x2).^2 - 1.5*x1 + 2.5*x2 + 1;
    
    end
     
end