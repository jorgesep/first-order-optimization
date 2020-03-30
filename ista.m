function [ fx_vals, x_100, x_next, resp ] = ista( A, b, x_true, x_0, opts )
%Implementation of ISTA algorithm


% set default parameters in case no given arguments.
if nargin == 0
    [A, b, x_true, x_0, opts] = default_parameters();
end

% parse input options, get default values in case not configured
[ t_k, lambda, maxiter, tol, verbose ] = parse_input_parameters(opts);

% console message.
fprintf('=====================================================\n')
fprintf(' ISTA: Iterative Shrinkage-Thresholding Algorithm\n')
fprintf('=====================================================\n')
fprintf(' Algorithm parameters:\n')
fprintf([' lambda: %1.3e ',...
         ' step  : %1.3e\n'],...
          lambda, t_k);

% get matrices
A_transpose = A';
A_square = A' * A;

% compute ||b|| only once
normb = norm(b); 

% initialize arrays 
ista_errors = zeros(4,maxiter);    % save some error as a stop criteria 
fx_vals = zeros(1,maxiter);        % objective function values

% initialize
x_k = x_0; niter = 0; ndiff = inf; 

% init function value
fx_k = lasso_function(A,x_k,b,lambda);

% true value of the objective function
fx_star = lasso_function(A,x_true,b,lambda); 

tic; 
% main loop
while and(ndiff>=tol, niter < maxiter)

    % gradient step
    x_grad = x_k - 2 * t_k .* (A_square * x_k - A_transpose* b);
    
    % soft thresholding
    x_next = subplus(abs(x_grad) - lambda * t_k) .* sign(x_grad);
    
    % function value for current iteration
    fx_next = lasso_function(A,x_next,b,lambda);
    
    % norm of the residual using x_k
    normr = norm(b - A*x_k);           % Save for printing
    if normr/normb < tol, break; end
    
    % error difference with previous values 
    ista_errors(1,niter+1) = norm(x_next - x_k)     /norm(x_k);
    ista_errors(2,niter+1) = norm(x_next - x_k,1)   /numel(x_k);
    ista_errors(3,niter+1) = norm(x_next - x_true,1)/norm(x_true,1);
    ista_errors(4,niter+1) = normr/normb ; 
    
    % display by console every 50 iterations
    if and(mod(niter,10) == 0, verbose )
        fprintf(['%d fx=[%.3e:%.3e:%3e:%3e] error=[%.3e:%.3e:%.3e] ',...
                 'residual=[%1.3e:%1.3e:%11.3e]\n'], ...
            niter, fx_next, (fx_k-fx_next), ...
            (fx_next-fx_star), (fx_next-fx_star)/fx_star, ...
            ista_errors(1,niter+1),ista_errors(2,niter+1), ...
            ista_errors(3,niter+1), normr, normr/normb,tol );
    end
    
    % update function vals 
    fx_k = fx_next;
    fx_vals(niter+1) = fx_k;
    ndiff = ista_errors(1,niter+1);

    % update
    x_k = x_next;
    niter = niter + 1;
    
    % Save vector of iteration 100
    if niter == 100 
        x_100 = x_next;
    end

end

resp = { 'ista' toc niter tol ndiff normr normb normr/normb };

save('ista_errors','ista_errors');

end

