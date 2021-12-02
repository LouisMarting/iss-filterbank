function ABCD_trxLine = trxLine_ABCD(Z0,g_eff,L)
    
    nF = length(g_eff);
    
    A = cosh(g_eff*L);
    B = Z0*sinh(g_eff*L); 
    C = Z0.^-1*sinh(g_eff*L);
    D = cosh(g_eff*L);
    
    A = reshape(A,[1,1,nF]);
    B = reshape(B,[1,1,nF]);
    C = reshape(C,[1,1,nF]);
    D = reshape(D,[1,1,nF]);
    
    AB           = cat(2,A,B);
    CD           = cat(2,C,D);
    ABCD_trxLine = cat(1,AB,CD);
    
    clear A B C D AB CD;
    
end