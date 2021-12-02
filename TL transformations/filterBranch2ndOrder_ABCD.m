function ABCD_eachFilterBranch = filterBranch2ndOrder_ABCD(Z_coup1,Z0_res,g_eff_res,L_res,Z_coup2,Z_coup_inter1)
    
    
    ABCD_coup1    = seriesLoad_ABCD(Z_coup1);
    ABCD_res1      = trxLine_ABCD(Z0_res,g_eff_res,L_res.*1.008);
    ABCD_res2      = trxLine_ABCD(Z0_res,g_eff_res,L_res.*1);
    ABCD_coup2    = seriesLoad_ABCD(Z_coup2); 
    ABCD_coup_inter1 = seriesLoad_ABCD(Z_coup_inter1); 
    
    ABCD_res_coup1        = mmat(ABCD_res1,ABCD_coup1,[1,2]);
    ABCD_res_1stOrder     = mmat(ABCD_coup_inter1,ABCD_res_coup1,[1,2]);
    ABCD_res_2filters     = mmat(ABCD_res2,ABCD_res_1stOrder,[1,2]);
    ABCD_eachFilterBranch = mmat(ABCD_coup2,ABCD_res_2filters,[1,2]);

end