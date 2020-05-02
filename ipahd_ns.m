function [ fx_k, x_100, x_next, resp ] = ipahd_ns( grad_fx, cost_fx, x_0, opts )
% IPAHD-NS: Inertial Proximal Algorithm with Hessian Damping Non-Smooth

% parse input options, get default values in case not configured
[ h, lambda, maxiter, tol, verbose ] = parse_input_parameters(opts);

% Global parameters of algorithm 
step = h ^ 2;

% damping coefficients (initial values taken as example from page 3.)
alpha = 3.1; % viscous 
%beta = 2*sqrt(step)*0.9; % Hessian-driven : b<2*sqrt(s)
beta = 0.02; % Hessian-driven : b<2*sqrt(s)
b_k = 1.0;

fprintf('=====================================================\n')
fprintf(' IPAHD-NS: Inertial Proximal Algorithm \n')
fprintf('           with Hessian Damping Non-Smooth\n')
fprintf('=====================================================\n')
fprintf(' lambda:%g\n step:%g\n alpha:%g\n beta:%g\n b:%g\n',... 
          lambda, step, alpha, beta, b_k);

% initialize arrays 
fx_k = zeros(1,maxiter); % objective function values
errors = zeros(4,maxiter);% array to keep criteria errors
x1x2 = zeros(2,maxiter) ; 
ix1 = opts.ix1; ix2 = opts.ix2 ;

% initialize
x_k = x_0; x_prev = x_0; k = 1; ndiff = inf; x_100 = x_0;

% load('x_true.mat'); load('A.mat'); load('b.mat');
% f0 = lasso_function(A,x_true,b,lambda);
% c0 = lasso_function(A,x_k,b,lambda);
% c1 = cost_fx(x_k,lambda);

tic;
% ipahd-ns main loop
while and(ndiff>=tol, k <= maxiter)
  
    % mu_k setup
    num_k = lambda * (k+alpha);
    den_k = lambda * (k+alpha) + k*(beta * sqrt(step) + step * b_k);
    mu_k = num_k/den_k ;

    % acceleration part
    alpha_k = 1 - alpha/(k+alpha);
    acc_k = alpha_k .* (x_k - x_prev);
    
    % gradient
    prox_lambda_f = prox_l1( x_k - grad_fx(x_k) , lambda) ; 
    Gfx_k = beta * sqrt(step)   * alpha_k * (1/lambda) * ( x_k - prox_lambda_f);
  
    % intermediate step 
    y_k = x_k + acc_k + Gfx_k ;

    % final step
    x_next = mu_k*y_k + (1-mu_k)*prox_l1(y_k-grad_fx(y_k),lambda/mu_k) ;

%     str1=sprintf('mu_k=%f alpha_k=%f acc_k=%f prox_lf=%f Gfx=%f y_k=%f x_next=%f cost=%f\n', ...
%         mu_k, alpha_k, norm(acc_k), norm(prox_lambda_f), norm(Gfx_k),norm(y_k),norm(x_next),cost_fx(x_next,lambda) );
%     fprintf('%s',str1);
    %--------------------------------------------------
    % end algorithm update
    %--------------------------------------------------
    
    % function value for current iteration
    fx_k(k) = cost_fx(x_next,lambda);
    
    % get criteria error for stoppping algorithm 
    [ ndiff, errors([1 2 3 4], k) ] = stop_criteria_value( 'norm2', x_next, x_k );

    % update
    x_prev = x_k; x_k = x_next;
 
    % save optimization value at iteration 100
    if k == 100, x_100 = x_next; end
    
    % save evolution of two variables (just for debugging)
    x1x2([1 2],k) = x_next([ix1 ix2]);
    
    % display messages every 10 iterations.
    if and(mod(k-1,10) == 0, verbose  )
        fprintf('%d cost = %.5s error = %1.6e\n', k, fx_k(k), ndiff);
    end
    k = k + 1;
    % show progress
    progressbar( maxiter, k, verbose )    

end
fprintf('\n');



% save most of configuration parameters in a cell
resp = { 'ipahd_ns' toc k lambda step tol ndiff  errors(2,end) };

% save useful files
nameStr = mfilename;
if ~exist(nameStr, 'dir'), mkdir(nameStr); end
save(strcat(nameStr,'/errors'),'errors');
save(strcat(nameStr,'/x1x2'),'x1x2');
save(strcat(nameStr,'/fx_k'),'fx_k');

writetable(struct2table(opts), strcat(nameStr,'/opts.txt'));


end
