function C = findCcoupler_Ql(series_or_parallel,resonatorLength,Ql,Qi,f0,Z1,Z2)
    
    R1 = real(Z1);
    R2 = real(Z2);
    X1 = imag(Z1);
    X2 = imag(Z2);
    
    if ~isfinite(Qi)
        Qc = 2.*Ql;
    else
        Qc = 2.*Qi.*Ql./(Qi-Ql);
    end
    
    if strcmpi('L4',resonatorLength)
        n_cycles = 2;
    elseif strcmpi('L2',resonatorLength)
        n_cycles = 1;
    end
    
    if strcmpi('SERIES',series_or_parallel)
        A = 1;
        B = 2 .* X1 + 2 .* X2;
        C = X1.^2 + 2 .* X1 .* X2 + X2.^2 - n_cycles .* Qc ./ pi .* 2 .* R1 .* R2 + (R1 + R2).^2;    
    elseif strcmpi('PARALLEL',series_or_parallel)
        A = (R1 + R2).^2 + (X1 + X2).^2 - n_cycles .* Qc ./ pi .* 2 .* R1 .* R2;
        B = 2 .* (R1 .* X2 + R2 .* X1) .* (R1 + R2) + 2 .*(X1 .* X2 - R1 .* R2) .* (X1 + X2);
        C = (R1 .* X2 + R2 .* X1).^2 + (X1 .* X2 - R1 .* R2).^2;
    end
    
    X = (-B - sqrt(B.^2 - 4 .* A .* C)) ./ (2 .* A);
    
    C = -1 ./ (2 .* pi .* f0 .* X);
end

