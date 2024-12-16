function fu = fisher_KPP_harvesting_nonlin(U, par)
    
    fu = par.r * U .* (1 - U) - par.b .* U;

end