function [ A, b, x, x_0, opts ] = default_parameters( )
%DEFAULT_ALGORITHM_PARAMETERS Initialize with default parameters

% load default image
if exist('cameraman.jpg','file') == 2
    img = im2double(imread('cameraman.jpg'));
    img = img(35:75,90:140);
else
    img = ...
   [0.7137    0.7098    0.6980    0.7176    0.7137    0.7098    0.7098    0.7137    0.7255    0.7137    0.7176;
    0.7059    0.7020    0.7098    0.7176    0.6941    0.7098    0.7098    0.7137    0.7216    0.7216    0.7255;
    0.7137    0.7020    0.7216    0.7176    0.7020    0.7176    0.7137    0.7059    0.7176    0.7216    0.7333;
    0.7137    0.7216    0.7216    0.7255    0.7176    0.7216    0.7176    0.7294    0.7176    0.7137    0.7098;
    0.7216    0.7255    0.7137    0.7333    0.7216    0.7255    0.7255    0.7294    0.7255    0.7216    0.7490;
    0.7176    0.7294    0.7216    0.7176    0.7294    0.7255    0.7333    0.7373    0.7333    0.7333    0.7333;
    0.7373    0.7098    0.7255    0.7098    0.7373    0.7373    0.7255    0.7333    0.7412    0.7294    0.7451;
    0.7059    0.7059    0.7098    0.7059    0.7098    0.7137    0.7176    0.7294    0.7333    0.7137    0.7098;
    0.7098    0.7255    0.7098    0.7176    0.7020    0.7098    0.7176    0.7294    0.7216    0.7333    0.6745;
    0.6902    0.7098    0.7059    0.7176    0.7059    0.7059    0.7216    0.7333    0.7294    0.7451    0.7333;
    0.7216    0.7216    0.7294    0.7176    0.7216    0.7294    0.7216    0.7412    0.7412    0.7529    0.7529];
end

% vectorize true image
x = img(:);

% blurring window
h = gaussian_kernel(9,4);
A = blur_operator(h,img);
L = max(eigs(A' * A)) ;
b = A * x; % blurring the image

% init img
x_0 = ones(length(b),1);

% general algorithm parameters 
% opts.step = 1 / L; 
% opts.lambda = 2e-5; 
% opts.maxiter = 100; 
% opts.tol = 1e-6;
% opts.verbose = false;

opts = generate_default_options;

end

