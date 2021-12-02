function L = findLcoupler(series_or_parallel,resonatorLength,Qt,Qi,f0,Z0_1,Z0_2)
    
    R0_1 = real(Z0_1);
    R0_2 = real(Z0_2);
    X0_1 = imag(Z0_1);
    X0_2 = imag(Z0_2);
    
    Qc = 2*(Qt.^-1 - Qi.^-1).^-1;  % Two couplers attached to the resonator
    
    w_0 = 2*pi*f0;
    
    if strcmpi('L4',resonatorLength)
        couplerEncountersPerCycle = 2;
    elseif strcmpi('L2',resonatorLength)
        couplerEncountersPerCycle = 1;
    end
    
    if strcmpi('SERIES',series_or_parallel)
        A = 1;
        B = 2*X0_1 + 2*X0_2;
        C = X0_1.^2 + 2*X0_1*X0_2 + X0_2.^2 - (couplerEncountersPerCycle.*Qc)./pi.*2.*R0_1.*R0_2 + (R0_1+R0_2).^2;    
    elseif strcmpi('PARALLEL',series_or_parallel)
        A = (R0_1+R0_2).^2 + (X0_1+X0_2).^2 - (couplerEncountersPerCycle.*Qc./pi.*2.*R0_1.*R0_2);
        B = 2.*(R0_1.*X0_2+R0_2.*X0_1).*(R0_1+R0_2) + 2.*(X0_1.*X0_2-R0_1.*R0_2).*(X0_1+X0_2);
        C = (R0_1.*X0_2+R0_2.*X0_1).^2 + (X0_1.*X0_2-R0_1.*R0_2).^2;
    end
    
    X = (-B + sqrt(B.^2 - 4.*A.*C))./(2.*A);
    
    L = X./w_0;
    
end