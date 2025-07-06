function Eta = Eta_q(etasd, Ebth, Tsd)
    Eta = etasd./DeltaH(sqrt(Ebth./Tsd));
end

function Delta = DeltaH(x)
    Delta = erf(x) - (2/sqrt(pi)).*x.*exp(-x.^2);
end