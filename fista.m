function [ fx_vals, x_100, x_next, resp ] = fista( A, b, x_true, x_0, opts )
% Particulat implementation of FISTA algorithm

% set default parameters in case no given arguments.
if nargin == 0
    [A, b, x_true, x_0, opts] = default_parameters();
end

% parse input options, get default values in case not configured
[ h, lambda, maxiter, tol, verbose ] = parse_input_parameters(opts);

% print out console message.
fprintf('=====================================================\n')
fprintf(' FISTA: Fast Iterative Shrinkage-Thresholding Algorithm\n')
fprintf('=====================================================\n')
fprintf(' Algorithm parameters:\n')
fprintf([' lambda: %1.3e\n',...
         ' step  : %1.3e\n'],...
          lambda, h);
%fprintf('=====================================================\n')

% get matrices
A_transpose = A';
A_square = A' * A;


% compute ||b|| only once
normb = norm(b);

% initialize arrays 
fista_errors = zeros(4,maxiter);    % criteria errors
fx_vals = zeros(1,maxiter); % objective function values

% initialize algorithm variables
x_k = x_0; y_k = x_0; t_k = 1; niter = 0; ndiff = inf;

% init function value
fx_k = lasso_function(A,x_k,b,lambda);

% true value of the objective function
value = lasso_function(A,x_true,b,lambda);

tic;
% main loop
while and(ndiff>=tol, niter < maxiter)
    
    % gradient step
    x_grad = y_k - 2 * h .* (A_square * y_k - A_transpose* b);
    
    % soft thresholding
    x_next = subplus(abs(x_grad) - lambda * h) .* sign(x_grad);

    % intermediate t step
    t_next = 0.5 * ( 1 + sqrt(1 + 4*t_k^2) );
    
    % FISTA update
    y_next = x_next + ((t_k - 1)/t_next) * (x_next - x_k);

    % function value for current iteration
    cost = lasso_function(A,y_next,b,lambda);
   
    residual = b - A*y_k;             % Compute residual with previous x
    normr = norm(residual);           % Save for printing
    if normr/normb < tol, break; end  % Test before solving; exit loop if true
        
    % error for stoppping criteria 
    fista_errors(1,niter+1) = norm(y_next - y_k)     /norm(y_k);
    fista_errors(2,niter+1) = norm(y_next - y_k,1)   /numel(y_k);
    fista_errors(3,niter+1) = norm(y_next - x_true,1)/norm(x_true,1);
    fista_errors(4,niter+1) = normr/normb ; 
    
    % display messages every 10 iterations.
    if and(mod(niter,10) == 0, verbose  )
        fprintf('%d cost = %.5s error = %1.6e\n', niter, cost, fista_errors(1,niter+1));
%         fprintf(['%d fx=[%.3e:%.3e:%3e:%3e] error=[%.3e:%.3e:%.3e] ', ...
%             'residual=[%1.3e:%1.3e:%1.3e]\n'], ...
%             niter, fx_next, (fx_k-fx_next), ...
%             (fx_next-value), (fx_next-value)/value, ...
%             fista_errors(1,niter+1),fista_errors(2,niter+1), ...
%             fista_errors(3,niter+1), ...
%             normr, normr/normb,tol);
    end
    
    % update function vals 
    fx_k = cost;
    fx_vals(niter+1) = fx_k;
    ndiff = fista_errors(1,niter+1);    
    
    % update
    x_k = x_next ;
    y_k = y_next ;
    t_k = t_next ;
    niter = niter + 1;    
    
    
    % save optimization value at iteration 100
    if niter == 100 
        x_100 = x_next;
    end

end

% save on a struct all main algorithm parameters
% resp.name = 'fista';
% resp.elapsed = toc;
% resp.niter = niter; resp.tol = tol; resp.ndiff = ndiff;
% resp.normr = normr; resp.normb = normb; resp.norme = normr/normb;
resp = { 'fista' toc niter tol ndiff normr normb normr/normb };

save('fista_errors','fista_errors');

end

