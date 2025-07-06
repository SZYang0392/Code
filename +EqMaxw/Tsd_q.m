function Tsd = Tsd_q(Ebth, Tsd, Etasd, Etasd_q)
    Tsd = 1.5*Etasd.*Tsd./Etasd_q./ThetaH(sqrt(Ebth./Tsd));
end

function Theta = ThetaH(x)
    Theta = 1.5*DeltaH(x) - (2./sqrt(pi))*x.^3.*exp(-x.^2);
end

function Delta = DeltaH(x)
    Delta = erf(x) - (2/sqrt(pi)).*x.*exp(-x.^2);
end