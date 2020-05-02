function [ opts ] = generate_default_options(  )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

L = 0.5 ;

% general algorithm parameters 
opts.step = 1 / L; 
opts.lambda = 2e-5; 
opts.maxiter = 100; 
opts.tol = 1e-6;
opts.verbose = false;

% parameters for debugging not needed for running algorithm
opts.ix1 = 1;
opts.ix2 = 2;

end

