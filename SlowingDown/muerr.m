function [mu, mup, mupp] = muerr(x)
    mu = gammainc(x, 1.5);
    mup = 2/sqrt(pi)*sqrt(x).*exp(-x);
    mupp = 2/sqrt(pi)*(1/2./sqrt(x) - sqrt(x)).*exp(-x);
end