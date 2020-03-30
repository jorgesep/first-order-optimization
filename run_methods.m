clear variables
close all

% algorithms to be tested
names = {'ista'; 'fista'} ;

% file names.
img_name = 'cameraman';
filename = strcat(img_name, '.jpg');
mat_file = strcat('A_',img_name, '.mat');
eig_file = strcat('A_eigens_', img_name, '.mat');

% init array struct with response parameters of run methods
header = {'name' 'elapsed' 'niter' 'tol' 'ndiff' 'normr' 'normb' 'norme'};
resp=cell(length(names),length(header));

% 2-D line plot of the optimization values in Y versus iterations.
display_plots = true;  % plot_results = 0; don't plot anything

% load image
img = im2double(imread(filename));

% vectorize true image
x = img(:);

% gaussian blurring kernel (gaussian)
h = gaussian_kernel(9,4); % window size nine and sigma 4.0

% Blurred linear operator
tic; fprintf('Getting linear operator ... ')
if exist('A.mat', 'file') == 2
    load('A.mat')
else
    A = blur_operator(h,img);
end
elapsed=toc; fprintf('elapsed %f\n', elapsed)

%*******************************************************************
% Local functions declaration.
%*******************************************************************
% calculate matrices
A_transpose = A';
A_square = A' * A;

% gradient
grad_fx = @(x) (A_square * x - A_transpose * b);

% blurring the image
tic; fprintf('Blurring image          ... ')
b = A * x;
elapsed=toc; fprintf('elapsed %f\n', elapsed)

% init with blurred image
x_0 = b;  

% maximum Eigen value
tic; fprintf('Computing eigenvalues   ... ')
if exist('Aeigs.mat', 'file') == 2
    load('Aeigs.mat')
else
    Aeigs = eigs(A' * A);
end
L = max(Aeigs);
elapsed = toc; fprintf('elapsed %f\n', elapsed)

% general algorithm parameters 
opts.step = 1 / L; 
opts.lambda = 2e-5; 
opts.maxiter = 3000; 
opts.tol = 1e-6;
opts.verbose = true;

% global configurations
opts.step = opts.step * 0.5;

% initialize array to save function outcome
costs  = zeros(length(names), opts.maxiter) ;

% run algorithms
for i=1:length(names)
switch names{i}
    case 'ista'

        [costs(i,:), x_mid, x_k, resp(i,:)]  = ista(A, b, x, x_0, opts);
        
    case 'fista'

        [costs(i,:), x_mid, x_k, resp(i,:)]  = fista(A, b, x, x_0, opts);

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
save('A','A');
save('Aeigs','Aeigs');

