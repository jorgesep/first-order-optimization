function [ prox_x ] = prox_l1( x, lambda )
% Soft-thresholding/shrinkage operator
% X = prox_l1(x,lambda) 


prox_x =  subplus(abs(x) - lambda) .* sign(x);

end

