function [ x_next, fobj ] = fista( A, b, x_0, step, lambda, maxiter, f_obj )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fprintf('=====================================================\n')
fprintf(' FISTA: Fast Iterative Shrinkage-Thresholding Algorithm\n')
fprintf('=====================================================\n')


% get matrices
tic; fprintf('Computing A_transpose * A ... ')
A_transpose = A';
A_square = A' * A;
elapsed = toc; fprintf('elapsed %f\n', elapsed)

% termination tolerance
tol = 1e-6;

% initialize
x_k = x_0; y_k = x_0; t_k = 1; niter = 0; ndiff = inf;


fobj = zeros(1,maxiter+1);
fobj(1) = f_obj(A,x_k,b,lambda);


while and(ndiff>=tol, niter < maxiter)
    % gradient step
    x_grad = y_k - 2 * step .* (A_square * y_k - A_transpose* b);
    
    % soft thresholding
    x_next = subplus(abs(x_grad) - lambda * step) .* sign(x_grad);
    
    % intermediate t step
    t_next = 0.5 * ( 1 + sqrt(1 + 4*t_k^2) );
    
    % FISTA update
    y_next = x_next + ((t_k - 1)/t_next) * (x_next - x_k);

    % error of stop criteria
    err_1 = norm(y_next - y_k, 1) / numel(y_next);
    err_2 = norm(y_next - y_k)/norm(y_k);
    t_step = (t_k - 1)/t_next ;
    sum_y_next = sum(y_next);
    
    % update
    x_k = x_next ;
    y_k = y_next ;
    t_k = t_next ;
    niter = niter + 1;

    % save function objective value
    fobj(niter+1) = f_obj(A,y_next,b,lambda);
    fdiff = fobj(niter) - fobj(niter+1);
    
    % display messages every 10 iterations.
    if mod(niter,10) == 0
    fprintf('%d fobj=%.6f fdiff=%.6e err_1=%3.5e err_2=%3.5e t_s=%3.3e %3.3f\n', ...
        niter, fobj(niter+1), fdiff, err_1, err_2, t_step, sum_y_next );
    end
    % save optimization value at iteration 100
    if niter == 100 
        x_100 = x_next;
    end

end

x_next = x_100;
end

