function [ fitness ] = testFunction( x, fhd, fNumber )
    
    global count;
    [~, counter] = size(x);
    count = count + counter;

    if strcmp(func2str(fhd), 'benchmarks')
        
        [dimension, pop] = size(x);
        fitness = zeros(1, pop);
        for i=1:pop
            fitness(i) = benchmarks(x(:,i), fNumber, dimension);
        end
        
    elseif strcmp(func2str(fhd), 'problems')
        
        [~, pop] = size(x);
        fitness = zeros(1, pop);
        for i=1:pop
            fitness(i) = problems(x(:,i), fNumber);
        end
        
    elseif strcmp(func2str(fhd), 'cec20_func_rw')
        
        fitness = feval(fhd, x', fNumber);
       
    else
        
        fitness = feval(fhd, x, fNumber);
                
    end

end
