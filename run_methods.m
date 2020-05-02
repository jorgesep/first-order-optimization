clear variables
close all

% algorithms to be tested
names = {'ipahd_ns'} ;
%names = {'fista'} ;
%names = {'ista'} ;
% names = {'igahd'} ;
names = {'ista'; 'fista'; 'igahd'; 'ipahd_ns'} ;


% file names.
img_name = 'cameraman';
filename = strcat(img_name, '.jpg');
mat_file = strcat('A_',img_name, '.mat');
eig_file = strcat('A_eigens_', img_name, '.mat');

% init array struct to save response of run methods
header = {'name' 'elapsed' 'niter' 'lambda' 'step' 'tol' 'error_2' 'error_1'};
resp=cell(length(names),length(header));

% 2-D line plot of the optimization values in Y versus iterations.
display_plots = false;  

% load image
img = im2double(imread(filename));

% vectorize true image
x = img(:);

% save true image in mat file
x_true=x; save('x_true','x_true'); clear x_true;  %#ok<NASGU>

% gaussian blurring kernel (gaussian)
h = gaussian_kernel(9,4); % window size nine and sigma 4.0

% Blurred linear operator
tic; fprintf('Getting linear operator ... ')
if exist('A.mat', 'file') == 2
    load('A.mat')
else
    A = blur_operator(h,img);
    save('A','A');
end
elapsed=toc; fprintf('elapsed %f\n', elapsed)

% calculate matrices
A_transpose = A';
A_square = A' * A;

% blurring the image
tic; fprintf('Blurring image          ... ')
b = A * x; save('b','b');
elapsed=toc; fprintf('elapsed %f\n', elapsed)

% init with blurred image
x_0 = b;  

% maximum Eigen value
tic; fprintf('Computing eigenvalues   ... ')
if exist('Aeigs.mat', 'file') == 2
    load('Aeigs.mat')
else
    Aeigs = eigs(A_square);
    save('Aeigs','Aeigs');
end
L = max(Aeigs);
elapsed = toc; fprintf('elapsed %f\n', elapsed)

%-------------------------------------------------------------------
% general algorithm parameters 
%-------------------------------------------------------------------
opts = generate_default_options; 
opts.step = 1e-2;
opts.lambda = 0.1;

opts.maxiter = 10000; 
opts.tol = 1e-6;
opts.verbose = false;

% global configurations

% Values presented in report number 1
opts.lambda = 2e-5;    % fista, ista
opts.step = 1 / (2*L); % fista, ista

opts.lambda = 0.001;

% two variables indexes to visualize algorithm evolution in a 2D plot.
opts.ix1 = 29630; opts.ix2=29379;
opts.ix1 = 32768; opts.ix2=59000;

% initialize array to save function outcome
costs  = zeros(length(names), opts.maxiter) ;

%-------------------------------------------------------------------
% Local functions declaration.
%-------------------------------------------------------------------

% gradient
%grad_fx = @(x) 2*(A_square * x - A_transpose * b);
grad_fx = @(x) 2 * A_transpose*(A*x - b);


% cost function
cost_fx = @(x,lambda) norm(A*x-b)^2 + lambda*sum(abs(x));

%-------------------------------------------------------------------
% main part: run algorithms
%-------------------------------------------------------------------
for i=1:length(names)
switch names{i}
    case 'ista'

        [costs(i,:), x_mid, x_k, resp(i,:)]  = ista(grad_fx, cost_fx, x_0, opts);
        
    case 'fista'

        [costs(i,:), x_mid, x_k, resp(i,:)]  = fista(grad_fx, cost_fx, x_0, opts);

    case 'igahd'

        [costs(i,:), x_mid, x_k, resp(i,:)]  = igahd(grad_fx, cost_fx, x_0, opts);
        
    case 'ipahd_ns'
        % parameters used in report 1
        opts.lambda = 0.001;  % ipahd_ns
        opts.step = 0.01;      % ipahd_ns
        opts.step = 0.02; % just testing
        [costs(i,:), x_mid, x_k, resp(i,:)]  = ipahd_ns(grad_fx, cost_fx, x_0, opts);

end
end


% Short scientific notation with 4 digits after the decimal point.
format shortE 

% display summary table 
disp(cell2table(resp, 'VariableNames',header))


if display_plots
    
    % plotting results
    plot_results(lasso_function(A,x,b,opts.lambda),names, 'log', costs);
        
end


% saving useful files



writetable(struct2table(opts), 'opts.txt')

