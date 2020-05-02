function progressbar( maxiter, idx, verbose )

unit_percent = floor(maxiter*0.02) ;


if and(mod(idx,unit_percent) == 0, ~verbose )
    if mod(floor(idx/(unit_percent)),5) == 0
        fprintf('|');
    else
        fprintf('.');
    end
end

end

