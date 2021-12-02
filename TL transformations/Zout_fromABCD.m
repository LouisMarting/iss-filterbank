function Zout = Zout_fromABCD(ABCD,Z_source)

    Z   = a2z(ABCD);
        
    Z11 = squeeze(Z(1,1,:));
    Z12 = squeeze(Z(1,2,:));
    Z21 = squeeze(Z(2,1,:));
    Z22 = squeeze(Z(2,2,:));
    
    Zout = Z22 - Z12.*Z21./(Z11+Z_source);
    Zout = Zout.';
    
end