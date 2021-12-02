function s = a2s(abcd,Z0)
    
    A = squeeze(abcd(1,1,:));
    B = squeeze(abcd(1,2,:));
    C = squeeze(abcd(2,1,:));
    D = squeeze(abcd(2,2,:));
    
    Z0_1 = Z0(:,1);
    Z0_2 = Z0(:,2);

    den = A.*Z0_2 + B + C.*Z0_1.*Z0_2 + D.*Z0_1;
    
    s(1,1,:) = (A.*Z0_2 + B - C.*conj(Z0_1).*Z0_2 - D.*conj(Z0_1))  ./den;
    s(1,2,:) = (2*(A.*D - B.*C).*(real(Z0_1).*real(Z0_2)).^0.5 )    ./den;
    s(2,1,:) = (2*(real(Z0_1).*real(Z0_2)).^0.5)                    ./den;
    s(2,2,:) = (-A.*conj(Z0_2) + B - C.*Z0_1.*conj(Z0_2) + D.*Z0_1) ./den;
    
end