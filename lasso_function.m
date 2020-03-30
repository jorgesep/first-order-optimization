function [ value ] = lasso_function( A,x,b,lambda )
% LASSO: Least Absolute Shrinkage and Selection Operator
% Local implementation of lasso function
% It computes difference of the blurred image with the A (linear operator)
% times a vector (true image) on the square 2-norm
% and the l1-norm regularization term applied on the true image
% regularized l1 norm of true values.

value = norm(b-A*x)^2 + lambda*sum(abs(x));

end

