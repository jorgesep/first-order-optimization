function [ e, errors ] = stop_criteria_value( type, y_next, y_k )
%return a value according norm2 or norm1 either.

    e_1 = norm(y_next-y_k)/norm(y_k);
    e_2 = norm(y_next-y_k,1)/numel(y_k);
    
    % verify A and final vector exists
    if exist('x_true','var') == 1
        e_3 = norm(y_next-x_true,1)/norm(x_true,1);
        e_4 = norm(y_next-x_true)/norm(x_true);
    else
        e_3 = inf;
        e_4 = inf;
    end

    switch type
        case 'norm2', e = e_1;
        case 'norm1', e = e_2;
        case 'xtrue', e = e_3;
        case 'residual', e = e_4;
        otherwise, e = inf;
    end
    
    errors = [e_1 e_2 e_3 e_4];

end

