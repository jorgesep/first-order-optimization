%*******************************************************************
% IPAHD: Inertial Proximal Algorithm with Hessian Damping
%*******************************************************************

clear all;
%*******************************************************************
% Local functions declaration.
%*******************************************************************

% objective function
f_obj = @(A, x, b, lambda) norm(A*x-b)^2 + lambda*sum(abs(x));

% Soft-thresholding/shrinkage operator
prox_l1 = @(x, lambda) subplus(abs(x) - lambda) .* sign(x);

%%
%*******************************************************************
% Initial configurations
%*******************************************************************
fprintf('=====================================================\n')
fprintf(' IPAHD-NS: Inertial Proximal Algorithm \n')
fprintf('           with Hessian Damping\n')
fprintf('=====================================================\n')


% name of files.
img_name = 'cameraman';
filename = strcat(img_name, '.jpg');
matfile = strcat('A_',img_name, '.mat');
eigfile = strcat('A_eigens_',img_name, '.mat');

% load image
img = im2double(imread(filename));

% vectorize true image
x = img(:);

% get blur kernel (gaussian) of window size 9 and sigma 1.5
h = gaussian_kernel(9,1.5);

% linear operator
tic; fprintf('Getting linear operator ... ')
if exist(matfile, 'file') == 2
    S = load(matfile,'A');
    A = S.A;
else
    A = blur_optimize(h,img);
    save(matfile,'A');
end
elapsed=toc; fprintf('elapsed %f\n', elapsed)

% transpose and square A matrix
tic; fprintf('Computing A_transpose times A ... ')
A_transpose = A';
A_square = A_transpose * A;
elapsed = toc; fprintf('elapsed %f\n', elapsed)

% step size (h) from maximum of eigen values
tic; fprintf('Computing eigenvalues ... ')

if exist(eigfile, 'file') == 2
    E = load(eigfile);
    A_square_eigs = E.A_square_eigs;
else
    A_square_eigs = eigs(A_square) ;
    save(eigfile,'A_square_eigs');
end
L = max(A_square_eigs);
step_size = 1 / L;
elapsed = toc; fprintf('elapsed %f\n', elapsed)

% blurring the image
tic; fprintf('Blurring image ... ')
b = A * x;
[brows, bcols] = size(b);
elapsed=toc; fprintf('elapsed %f\n', elapsed)

%%
%*******************************************************************
% IPAHD algorithm part.
%*******************************************************************


% Global parameters of algorithm 
lambda = 0.001;
step_size = 0.01;
step = step_size ^ 2;
number_iterations = 500;

% damping coefficients (initial values taken as example from page 3.)
alpha = 3.1; % viscous 
beta = 0.02; % Hessian-driven : b<2*sqrt(s)
b_k = 1;

fprintf('=====================================================\n')
fprintf(' Algorithm parameters:')
fprintf([' lambda: %f\n',...
         ' L     : %f\n',...
         ' step  : %f\n',...
         ' alpha : %f\n',...
         ' beta  : %f\n',...
         ' b     : %f\n'],... 
          lambda, L, step, alpha, beta, b_k);
fprintf('=====================================================\n')

% init with blurred image
x_k = b; 
x_prev = b;


% vector with values of the objective function
fval = zeros(1,number_iterations+1);
fval(1) = f_obj(A,x_k,b,lambda);
total_elapsed = 0;

% ipahd main loop
for k=1:number_iterations
    
    tic; 
    
    % k threshold argument for the proximal operator
    mu_k = k/(k+alpha) * (beta * sqrt(step) + step * b_k);
    
    % derivative wrt time
    Dfx_k = (1 - alpha/(k+alpha)) .* (x_k - x_prev);
    
    % gradient
    x_grad = x_k - 2 * (1 / L) .* (A_square * x_k - A_transpose* b);
    Gfx_k = beta * sqrt(step) * (1 - alpha/(k+alpha)) * x_grad ;
    
    % intermediate step (acceleration path)
    y_k = x_k + Dfx_k + Gfx_k ;
    
    % proximal operator: soft thresholding - shrinkage
    x_next = subplus(abs(y_k) - mu_k) .* sign(y_k) ;
 
    elased = toc;
    err_2 = norm(x_next - x_prev)/norm(x_prev);
    
    % updating
    x_prev = x_k;
    x_k = x_next;
    
    % evaluating the objective function
    fval(k+1) = f_obj(A,x_k,b,lambda);
    ftol = fval(k) - fval(k+1);
    fprintf('%d fobj=%f ftol=%.6f err_2= %.6f %3.3e\n', ...
    k, fval(k+1), ftol, err_2, elapsed )

  %if abs(ftol) < 1e-4
  %    break;
  %end
end

prImg = reshape(x_k,size(img));
figure(1); imshow(im2uint8(prImg));title('Recovered');


