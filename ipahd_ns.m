function [ fx_vals, x_100, x_next, resp ] = ipahd_ns( A, b, x_true, x_0, opts )
% IPAHD-NS: Inertial Proximal Algorithm with Hessian Damping Non-Smooth

%*******************************************************************
% Initial configurations
%*******************************************************************

% set default parameters in case no given arguments.
if nargin == 0
    [A, b, x_true, x_0, opts] = default_parameters();
end

% parse input options, get default values in case not configured
[ step, lambda, maxiter, tol, verbose ] = parse_input_parameters(opts);

fprintf('=====================================================\n')
fprintf(' IPAHD-NS: Inertial Proximal Algorithm \n')
fprintf('           with Hessian Damping Non-Smooth\n')
fprintf('=====================================================\n')

% get matrices
A_transpose = A';
A_square = A' * A;

% Global parameters of algorithm 
step_size = 0.01;
step = step_size ^ 2;

%*******************************************************************
% Local functions declaration.
%*******************************************************************

% gradient
grad_fx = @(x) 2 * step * (A_square * x - A_transpose* b);

%%
%*******************************************************************
% IPAHD-NS algorithm part.
%*******************************************************************

% damping coefficients (initial values taken as example from page 3.)
alpha = 3.1; % viscous 
%beta = 2*sqrt(step)*0.9; % Hessian-driven : b<2*sqrt(s)
beta = 0.02; % Hessian-driven : b<2*sqrt(s)
b_k = 1;

fprintf(' Algorithm parameters:\n')
fprintf([' lambda:%1.3e',...
         ' step:%1.3e',...
         ' alpha:%1.3e',...
         ' beta:%1.3e',...
         ' b:%f\n'],... 
          lambda, step, alpha, beta, b_k);

% init function value
fx_k = lasso_function(A,x_0,b,lambda);

% expected value of the function
value = lasso_function(A,x_true,b,lambda);
      
% compute ||b|| only once
normb = norm(b);

% initialize arrays 
ipahd_ns_errors = zeros(4,maxiter);    % criteria errors
fx_vals = zeros(1,maxiter); % objective function values

% initialize
x_k = x_0; x_prev = x_0; k = 0; ndiff = inf; 


%%
tic;
% ipahd-ns main loop
while and(ndiff>=tol, k <= maxiter)
 fprintf('0 tol:%1.3e\n', tol);   
    %--------------------------------------------------
    % init algorithm step
    %--------------------------------------------------
    
    % mu_k setup
    num_k = lambda * (k+alpha);
    den_k = lambda * (k+alpha) + k*(beta * sqrt(step) + step * b_k);
    mu_k  = num_k/den_k;
    
    % first derivative wrt time
    alpha_k = 1 - alpha/(k+alpha);
    Dfx_k = alpha_k .* (x_k - x_prev);
    
    % gradient
    prox_lambda = prox_l1( x_k - grad_fx(x_k) , lambda) ;
    Gfx_k = beta * sqrt(step) * (1/lambda) * alpha_k * ( x_k - prox_lambda);
    
    % intermediate step (acceleration path)
    y_k = x_k + Dfx_k + Gfx_k ;
    
    % final step
    x_next = mu_k*y_k + (1-mu_k)*prox_l1(y_k-grad_fx(y_k),lambda/mu_k) ;
    
    %--------------------------------------------------
    % end algorithm update
    %--------------------------------------------------
    
    % final value for current iteration
    fx_vals(k+1) = lasso_function(A,x_next,b,lambda);
    
    
    residual = b - A*x_k;             % Compute residual with previous x
    normr = norm(residual);           % Save for printing
    if normr/normb < tol, break; end  % Test before solving; exit loop if true
        
    
    
    % error for stoppping criteria 
    ipahd_ns_errors(1,k+1) = norm(x_next - x_k)     /norm(x_k);
    ipahd_ns_errors(2,k+1) = norm(x_next - x_k,1)   /numel(x_k);
    ipahd_ns_errors(3,k+1) = norm(x_next - x_true,1)/norm(x_true,1);
    ipahd_ns_errors(4,k+1) = normr/normb ; 
    ndiff = ipahd_ns_errors(1,k+1);

    % updating
    x_prev = x_k;
    x_k = x_next;
    
    
    fprintf(['%d fx=[%1.3e:%1.3e:%1.3e:%1.3e] ', ...
                'error=[%1.3e:%1.3e:%1.3e] ', ...
                'residual=[%1.3e:%1.3e:%1.3e:%1.3e]\n'], ...
             k, fx_vals(k+1), fx_k - fx_vals(k+1), fx_vals(k+1)-value, ...
             (fx_vals(k+1)-value)/value, ipahd_ns_errors(1,k+1), ipahd_ns_errors(2,k+1), ...
             ipahd_ns_errors(3,k+1),  normr, normr/normb, tol )
    
    % display messages every 10 iterations.
    %if and(mod(k,10) == 0, verbose  )
%        fprintf(['%d fx=[%1.3e:%1.3e:%1.3e:%1.3e] error=[%1.3e:%1.3e:%1.3e] residual=[%1.3e:%1.3e:%1.3e]\n'], ...
%            k, fx_vals(k), (fx_k - fx_vals(k)), (fx_vals(k)-value), (fx_vals(k)-value)/value, ipahd_ns_errors(1,k),ipahd_ns_errors(2,k), ...
%            ipahd_ns_errors(3,k),  normr, normr/normb,tol);
   % end
    
    k = k + 1;
    fx_k = fx_vals(k);
    
    % Save vector of iteration 100
    if k == 100 
        x_100 = x_next;
    end
fprintf('1 tol:%1.3e\n', tol);
end

% save on a struct all main algorithm parameters
resp.name = 'ipahd_ns';
resp.elapsed = toc;
resp.niter = k; resp.tol = tol; resp.ndiff = ndiff;
resp.normr = normr; resp.normb = normb; resp.norme = normr/normb;

save('ipahd_ns_errors','ipahd_ns_errors');


end
