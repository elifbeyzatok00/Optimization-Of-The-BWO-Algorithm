
    rand('state', sum(100*clock));
    if strcmp(func2str(fhd), 'benchmarks')
        [lowerBand, upperBand] = terminate_benchmark(fNumber);
        lbArray = ones(1, dimension) * lowerBand;
        ubArray = ones(1, dimension) * upperBand;
    elseif strcmp(func2str(fhd), 'problems')
        [lbArray, ubArray] = terminate_problem(fNumber);
    elseif strcmp(func2str(fhd), 'cec20_func_rw')
        [par] = Cal_par(fNumber);
        lbArray = par.xmin; ubArray = par.xmax;
    else
        lbArray = ones(1, dimension) * -100;
        ubArray = ones(1, dimension) * 100;
    end