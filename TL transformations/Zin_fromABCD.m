function Zin = Zin_fromABCD(ABCD,Z_load)

    Z   = a2z(ABCD);
        
    Z11 = squeeze(Z(1,1,:));
    Z12 = squeeze(Z(1,2,:));
    Z21 = squeeze(Z(2,1,:));
    Z22 = squeeze(Z(2,2,:));
    
    Zin = Z11 - Z12.*Z21./(Z22+Z_load);
    Zin = Zin.';
    
end