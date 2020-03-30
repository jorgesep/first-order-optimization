function [ prox_x ] = prox_l1( x, lambda )
% X = prox_l1(x,lambda) 
% Soft-thresholding/shrinkage operator


prox_x =  subplus(abs(x) - lambda) .* sign(x);

end

