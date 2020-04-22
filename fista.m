function [ fx_k, x_100, y_next, resp ] = fista( A, b, x_true, x_0, opts )
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
fprintf([' lambda: %1.3f\n',...
         ' step  : %1.3f\n'],...
          lambda, h);

% get matrices
A_transpose = A';
A_square = A' * A;

% gradient
grad_fx = @(x) (A_square * x - A_transpose * b);

% compute ||b|| only once
normb = norm(b);

% initialize arrays 
fx_k = zeros(1,maxiter); % objective function values
x1x2 = zeros(2, maxiter) ; ix1=29630; ix2=29379;

% initialize algorithm variables
x_k = x_0; y_k = x_0; t_k = 1; niter = 0; ndiff = inf;


tic;
% main loop
while and(ndiff>=tol, niter < maxiter)
    
    % gradient step
    x_grad = y_k - 2 * h .* grad_fx(y_k);
    
    % soft thresholding
    x_next = prox_l1(x_grad,lambda * h);

    % intermediate t step
    t_next = 0.5 * ( 1 + sqrt(1 + 4*t_k^2) );
    
    % FISTA update
    y_next = x_next + ((t_k - 1)/t_next) * (x_next - x_k);

    % update counter
    niter = niter + 1;
    
    % function value for current iteration
    fx_k(niter) = lasso_function(A,y_next,b,lambda);

    % error for stoppping criteria 
    ndiff = get_stop_criteria_val('norm2');

    % update
    x_k = x_next ; y_k = y_next ; t_k = t_next ;

    % save optimization value at iteration 100
    if niter == 100, x_100 = y_next; end
    
    % display messages every 10 iterations.
    if and(mod(niter-1,10) == 0, verbose  )
        fprintf('%d cost = %.5s error = %1.6e\n', niter, fx_k(niter), errors(1,niter));
    end
    % show progress
    progressbar( maxiter, niter, verbose )
    
end
fprintf('\n');


    
% save all in a cell
resp = { 'fista' toc niter lambda h tol ndiff  errors(2,end) };

nameStr = mfilename;
if ~exist(nameStr, 'dir'), mkdir(nameStr); end
save(strcat(nameStr,'/errors'),'errors');
save(strcat(nameStr,'/x1x2'),'x1x2');

% local function to save in an array different versions error stop criteria
function e = get_stop_criteria_val(type)
    % create array
    if niter == 1 
        errors = zeros(4,maxiter);    % array for keep criteria errors
    end
    e_1 = norm(y_next-y_k)/norm(y_k);
    e_2 = norm(y_next-y_k,1)/numel(y_k);
    e_3 = norm(y_next-x_true,1)/norm(x_true,1);
    e_4  = norm(b - A*y_k)/normb;
    errors([1 2 3 4],niter) = [e_1 e_2 e_3 e_4];

    x1x2([1 2],niter) = y_next([ix1 ix2]);
    switch type
        case 'norm2'
            e = errors(1,niter);
        case 'norm1'
            e = errors(2,niter);
        case 'xtrue'
            e = errors(3,niter);
        case 'residual'
            e = errors(4,niter);
        otherwise
            e = inf;
    end

end
end

