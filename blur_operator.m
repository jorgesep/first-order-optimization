function [ A ] = blur_operator( h, img )
%BLUR_OPTIMIZE Summary of this function goes here
%   Detailed explanation goes here

%

[hrows, hcols] = size(h);
[irows, icols] = size(img);

h_size =  hrows * hcols;
img_size = irows * icols; 

% number of nonzero elements in matrix A.
nnz_elements = img_size * h_size; 

% vectors to build sparse matrix
row_indexes = 1:img_size;
row_indexes = repelem(row_indexes, h_size);
col_indexes = zeros(1, nnz_elements);
values      = zeros(1, nnz_elements);

% create a zero matrix with the same size of input image
h_img = zeros(irows, icols);

% copy window kernel in position (1,1) of zero matrix.
h_img(1:hrows,1:hcols) = h;

% placing center of h window in the first row/colum of zero matrix 
h_img = circshift(h_img,[-floor(hrows/2),-floor(hcols/2)]);

% shift h window
for j=0:icols-1

    col_shifted = circshift(h_img, [0, j]);

    for i=0:irows-1
        h_shifted = circshift(col_shifted, [i,0]);

   
        % get indexes and values
        [~,col,v] = find(h_shifted(:)');

        h_index = (irows * j) + i;

        idx = (h_index * h_size) + 1 :  (h_index * h_size) + h_size;

        %fprintf('j=%d i=%d vector_index=%d\n',j,i,h_index) 
        col_indexes( idx ) = col;
        values( idx ) = v;
    end
end

A = sparse(row_indexes',col_indexes',values');


end

