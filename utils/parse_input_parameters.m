function [ step, lambda, maxiter, tol, verbose ] = parse_input_parameters( opts )
%Auxiliar function to parse the algorithm input parameters.

if ~isfield(opts, 'maxiter')
    maxiter = 100;
else
    maxiter = opts.maxiter;
end

if ~isfield(opts, 'step')
    step = 0.1 ;
else
    step = opts.step;
end

if ~isfield(opts, 'lambda')
    lambda = 2e-5 ;
else
    lambda = opts.lambda;
end

if ~isfield(opts, 'tol')
    tol = 1e-6 ;
else
    tol = opts.tol;
end

if ~isfield(opts, 'verbose')
    verbose = false ;
else
    verbose = opts.verbose;
end


end

