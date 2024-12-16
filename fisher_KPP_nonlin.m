function fu = fisher_KPP_nonlin(U, par)
    
    fu = par.r .* U .* (1-U);

end