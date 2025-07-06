function tau = Tause(Ps, n)
% Calculate fast ion slowing-down time by electrons
    epsilon_0 = 8.854187817e-12;
    e = 1.602176565e-19;

    lnL = LnLambda(Ps, 1, n);
    Gammane = (e^4/4/pi/(epsilon_0^2))*(Ps(n).q/Ps(n).m)^2*Ps(1).n0.*lnL;
    vte = sqrt(2*Ps(1).T*e/Ps(1).m);
    tau = 3*sqrt(pi)/4*Ps(1).m/Ps(n).m./Gammane.*(vte.^3);
end