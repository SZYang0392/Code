function lnL = LnCoulomb(Ps, Z)
% Calculate Coulomb logarithm approximately
    lnL = 30.4 + 1.5*log(Ps(1).T) - 0.5*log(Ps(1).n0) - log(Z);
end