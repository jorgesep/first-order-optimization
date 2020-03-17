%*******************************************************************
% IGAHD-NS: Inertial Gradient Algorithm with Hessian Damping
%*******************************************************************


clear;


%%
%*******************************************************************
% Initial configurations
%*******************************************************************
fprintf('=====================================================\n')
fprintf(' IGAHD: Inertial Gradient Algorithm \n')
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
h = gaussian_kernel(9,4);

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
% Local functions declaration.
%*******************************************************************

% objective function
f_obj = @(A, x, b, lambda) norm(A*x-b)^2 + lambda*sum(abs(x));


% Soft-thresholding/shrinkage operator
prox_l1 = @(x, lambda) subplus(abs(x) - lambda) .* sign(x);

% gradient
Grad_fx = @(x) 2 * (A_square * x - A_transpose* b);

%%
%*******************************************************************
% IGAHD algorithm part.
%*******************************************************************


% Global parameters of algorithm 
lambda = 0.001;
step_size = 0.01;
step = step_size ^ 2;
%step = 1 / L;
number_iterations = 1000;

% damping coefficients (initial values taken as example from page 3.)
alpha = 3.1; % viscous 
beta = 0.02; % Hessian-driven : b<2*sqrt(s)
%beta = 0;

fprintf('=====================================================\n')
fprintf(' Algorithm parameters:')
fprintf([' lambda: %f\n',...
         ' L     : %f\n',...
         ' step  : %f\n',...
         ' alpha : %f\n'],...
          lambda, L, step, alpha );
fprintf('=====================================================\n')


% init with blurred image
x_k = b; 
x_prev = b;


% vector with values of the objective function
fval = zeros(1,number_iterations+1);
fval(1) = f_obj(A,x_k,b,lambda);

% f_obj(x_min): function evaluated with the true image
fx_min = f_obj(A,x,b,lambda);
total_elapsed = 0;

% main loop of the ipahd algorithm
for k=1:number_iterations
    
    tic; 
    
    alpha_k = 1 - alpha/k;
    
    % derivative part
    Dfx_k = alpha_k .* (x_k - x_prev);
    
    % Hessian part. See P.10 
    Hfx_k = beta * sqrt(step)* (Grad_fx(x_k) - Grad_fx(x_prev));
    
    % gradient
    Gfx_k = beta * sqrt(step) * (1/k) * Grad_fx(x_prev) ;
    
    y_k = x_k + Dfx_k - Hfx_k - Gfx_k;
    
    % final step
    x_next = y_k - step * Grad_fx(y_k);

    elapsed = toc;
    total_elapsed = total_elapsed + elapsed ;
    err_2 = norm(x_next - x_prev)/norm(x_prev);
    err_1 = norm(x_next - x_prev,1)/numel(x_next);
    
    % updating
    x_prev = x_k;
    x_k = x_next;
    
    % comparing values of the objective function
    fval(k+1) = f_obj(A,x_k,b,lambda);
    ftol = fval(k) - fval(k+1);
    
    if mod(k,10) == 0
    fprintf('%d fobj=%f ftol=%.9f err_1=%3.3e err_2=%3.3e %3.3e\n', ...
    k, fval(k+1), ftol, err_1, err_2, elapsed )    
    end
    % Save vector of iteration 100
    if k == 100 
        x_100 = x_next;
    end


  %if abs(ftol) < 1e-4
  %    break;
  %end
end

prImg = reshape(x_k,size(img));
figure(1); imshow(im2uint8(prImg));title('Recovered');

fx_min = f_obj(A,x,b,lambda);
figure(2);
vals = fval-fx_min ;
results = vals(vals > 0) ;

loglog(results);
ylabel('${f(x_k)-f(x^{*})}$','interpreter','latex', 'FontWeight','bold')
xlabel('${k}$','interpreter','latex', 'FontWeight','bold')
title('${l_1}$','interpreter','latex', 'FontWeight','bold')
legend({'IGAHD'},'Location','northeast')
grid on;


