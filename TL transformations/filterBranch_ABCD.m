function ABCD_eachFilterBranch = filterBranch_ABCD(Z_coup1,Z0_res,g_eff_res,L_res,Z_coup2)
    
    
    ABCD_coup1    = seriesLoad_ABCD(Z_coup1);
    ABCD_res      = trxLine_ABCD(Z0_res,g_eff_res,L_res);
    ABCD_coup2    = seriesLoad_ABCD(Z_coup2);  
    
    ABCD_res_coup1        = mmat(ABCD_res,ABCD_coup1,[1,2]);
    ABCD_eachFilterBranch = mmat(ABCD_coup2,ABCD_res_coup1,[1,2]);

end