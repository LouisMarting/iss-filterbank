function z = a2z(abcd)
    
    A = abcd(1,1,:);
    B = abcd(1,2,:);
    C = abcd(2,1,:);
    D = abcd(2,2,:);
    
    z(1,1,:) = A./C;
    z(1,2,:) = (A.*D-B.*C)./C;
    z(2,1,:) = 1./C;
    z(2,2,:) = D./C;
    
end