function Z = ztransfer(ZL,Z0,k,l)
    Z                = Z0 .* (ZL+1i*Z0.*tan(k*l)) ./ (Z0+1i*ZL.*tan(k*l));
    Z(~isfinite(ZL)) = -1i*Z0.*cot(k(~isfinite(ZL))*l);
    Z(isnan(k))      = Inf;
    Z(k==0)          = Inf;
end