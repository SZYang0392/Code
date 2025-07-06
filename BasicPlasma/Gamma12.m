function G = Gamma12(P1, P2)
    e = 1.602176565e-19;
    epsilon_0 = 8.854187817e-12;
    G = P2.n0.*(P1.q.*P2.q.*e^2).^2./(4*pi*epsilon_0^2*P1.m^2).*LnLambda([P1, P2], 1, 2);
end