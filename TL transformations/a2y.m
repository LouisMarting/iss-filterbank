function y = a2y(abcd)
    
    A = abcd(1,1,:);
    B = abcd(1,2,:);
    C = abcd(2,1,:);
    D = abcd(2,2,:);
    
    y(1,1,:) = D./B;
    y(1,2,:) = (B.*C-A.*D)./B;
    y(2,1,:) = -1./B;
    y(2,2,:) = A./B;
end