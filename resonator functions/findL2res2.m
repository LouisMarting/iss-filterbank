function [lres]=findL2res2(Z_0,Z_KID,Z_0_f,ep,f0,C_coup1,C_coup2)
    % Find lres forcing Z0/2 as input impedance of the resonator


    %%
    mu_0   = pi*4e-7;             % []    Permeability of FS
    eps_0  = 8.854187817620e-12;  % []    Permittivity of FS
    c_0    = 1/sqrt(eps_0*mu_0);  % [m/s] Speed of light in FS


    %% Filter
    w0 = 2*pi*f0;
    v = c_0/sqrt(ep);
    k0 = w0/v;
    
    Z_coup1 = -1i./(w0.*C_coup1);

    Z_coup2 = -1i./(w0.*C_coup2);

    A = 0;
    
    kl = atan((Z_0/2-Z_coup2-A)./(-1i*(Z_0/2./Z_0_f.*A-Z_coup2./Z_0_f.*A-Z_0_f)));
    kl(kl<0) = kl(kl<0)+pi;
    
    lres = real(kl./k0);

end