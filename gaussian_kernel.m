function [ h ] = gaussian_kernel( arg1, arg2 )
%GAUSSIAN_KERNEL gaussian window
%   This generates and return a matrix with gaussian values
  
  % window size
  if (nargin > 0 && isreal(arg1))
    if mod(arg1,2) == 0
      hsize = arg1 + 1;
    else
      hsize = arg1 ;
    end
  else
    hsize = 3;
  end

  % Get sigma
  if (nargin > 1 && isreal(arg2) && length(arg2(:)) == 1)
    sigma = arg2;
  else
    sigma = 0.5;
  end
  
  % center to zero the array
  middle = floor(hsize/2);
  x = -middle:1:middle;
  
  % Compute Gaussian
  p1 = -.5 * ((x/sigma) .^ 2);
  p2 = (sigma * sqrt(2*pi)); 
  f = exp(p1) ./ p2;
  f = f ./ sum(f);
  h = f'*f;
end

