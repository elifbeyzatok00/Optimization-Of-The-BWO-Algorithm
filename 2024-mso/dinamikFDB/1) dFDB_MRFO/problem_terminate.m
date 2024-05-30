function [ecosize, n, maxFE, lb, ub] = problem_terminate()

    % Parameter settings:
    
    % ecosystem (population) size
    ecosize = 50;

    % problem dimension
    n = 30;

    % number of function evaluations
    maxFE = 10000 * n;

    % problem lower band 
    lowerBand = -100;
    lb = ones(1, n) * lowerBand;

    % problem upper band
    upperBand = 100;
    ub = ones(1, n) * upperBand;

end

