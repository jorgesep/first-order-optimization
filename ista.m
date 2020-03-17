function [ x_next, fobj ] = ista( A, b, x_0, step, lambda, maxiter, f_obj )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fprintf('=====================================================\n')
fprintf(' ISTA: Iterative Shrinkage-Thresholding Algorithm\n')
fprintf('=====================================================\n')


% get matrices
tic; fprintf('Computing A_transpose * A ... ')
A_transpose = A';
A_square = A' * A;
elapsed = toc; fprintf('elapsed %f\n', elapsed)

% termination tolerance
tol = 1e-6;

% initialize
x_k = x_0; niter = 0; ndiff = inf;


fobj = zeros(1,maxiter+1);
fobj(1) = f_obj(A,x_k,b,lambda);


while and(ndiff>=tol, niter < maxiter)
    % gradient step
    x_grad = x_k - 2 * step .* (A_square * x_k - A_transpose* b);
    
    % soft thresholding
    x_next = subplus(abs(x_grad) - lambda * step) .* sign(x_grad);
    
    % difference with previous values 
    ndiff = norm(x_next - x_k)/norm(x_k);
    
    % update
    x_k = x_next;
    niter = niter + 1;
    
    % save function objective value
    fobj(niter+1) = f_obj(A,x_next,b,lambda);
    fdiff = fobj(niter) - fobj(niter+1);
    
    % display by console every 10 iterations
    if mod(niter,10) == 0
    fprintf('%d fobj=%.6f fdiff=%.6e ndiff=%.6f\n', ...
        niter, fobj(niter+1), fdiff, ndiff );
    end
    
    % Save vector of iteration 100
    if niter == 100 
        x_100 = x_next;
    end
end

x_next = x_100;
end

