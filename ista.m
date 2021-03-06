function [ fx_k, x_100, x_next, resp ] = ista( grad_fx, cost_fx, x_0, opts )
%Implementation of ISTA algorithm


% parse input options, get default values in case not configured
[ step, lambda, maxiter, tol, verbose ] = parse_input_parameters(opts);

% console message.
fprintf('=====================================================\n')
fprintf(' ISTA: Iterative Shrinkage-Thresholding Algorithm\n')
fprintf('=====================================================\n')
fprintf(' Algorithm parameters:\n')
fprintf([' lambda: %g\n',...
         ' step  : %g\n'],...
          lambda, step);
      
% initialize arrays 
fx_k = zeros(1,maxiter);  % objective function values
errors = zeros(4,maxiter);% array to keep stop criteria errors
x1x2 = zeros(2,maxiter) ; 
ix1 = opts.ix1; ix2 = opts.ix2 ;

% initialize algorithm variables
x_k = x_0; niter = 0; ndiff = inf; x_100 = x_0;

tic; 
% main loop
while and(ndiff>=tol, niter < maxiter)

    % gradient step
    x_grad = x_k - step .* grad_fx(x_k);

    % soft thresholding
    x_next = prox_l1(x_grad,lambda * step);

    % update counter
    niter = niter + 1;
    
    % function value for current iteration
    fx_k(niter) = cost_fx(x_next,lambda);
    
    % get criteria error for stoppping algorithm 
    [ ndiff, errors([1 2 3 4], niter) ] = stop_criteria_value( 'norm2', x_next, x_k );
    
    % update
    x_k = x_next ;
    
    % save optimization value at iteration 100
    if niter == 100, x_100 = x_next; end
    
    % save evolution of two variables (just for debugging)
    x1x2([1 2],niter) = x_next([ix1 ix2]);
    
    % display messages every 10 iterations.
    if and(mod(niter-1,10) == 0, verbose  )
        fprintf('%d cost = %.5s error = %1.6e\n', niter, fx_k(niter), ndiff);
    end
    
    % show progress
    progressbar( maxiter, niter, verbose )
    
end
fprintf('\n');
    
% save most of configuration parameters in a cell
resp = { 'ista' toc niter lambda step tol ndiff  errors(2,end) };

% save useful files
nameStr = mfilename;
if ~exist(nameStr, 'dir'), mkdir(nameStr); end
save(strcat(nameStr,'/errors'),'errors');
save(strcat(nameStr,'/x1x2'),'x1x2');
save(strcat(nameStr,'/fx_k'),'fx_k');

writetable(struct2table(opts), strcat(nameStr,'/opts.txt'));

end

