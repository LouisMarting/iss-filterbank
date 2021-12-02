function ABCD = seriesLoad_ABCD(Z_load)
    
    nF = length(Z_load);
    
    A    = ones(1,nF);
    B    = Z_load;
    C    = zeros(1,nF);
    D    = ones(1,nF);
    
    A    = reshape(A,[1,1,nF]);
    B    = reshape(B,[1,1,nF]);
    C    = reshape(C,[1,1,nF]);
    D    = reshape(D,[1,1,nF]);
    
    AB   = cat(2,A,B);
    CD   = cat(2,C,D);
    ABCD = cat(1,AB,CD);

end