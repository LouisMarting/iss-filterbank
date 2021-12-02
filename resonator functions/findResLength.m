function [lres]=findResLength(f0,lam0,Z1,Z2,Zres,C1,C2)
    % Find lres forcing Z0/2 as input impedance of the resonator
    
    if Z1 == 0
        Z_C1 = 0;
    else
        Z_C1 = -1i ./ (2 .* pi .* f0 .* C1);
    end
    if Z2 == 0
        Z_C2 = 0;
    else
        Z_C2 = -1i ./ (2 .* pi .* f0 .* C2);
    end

    k0 = 2 .* pi ./ lam0;

    A = Z_C2 + Z2;
    
    kl = atan( (Z1 - Z_C1 - A) ./ (-1i .* (Z1 .* A ./ Zres - Z_C1 .* A ./ Zres - Zres)));
    kl(kl<0) = kl(kl<0) + pi;
    
    lres = real(kl ./ k0);
end