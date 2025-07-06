function lnL = LnLambda(Ps, m, n)
% Calculate Coulomb logarithm
    e = 1.602176565e-19;
    epsilon_0 = 8.854187817e-12;

    Ld = Debye(Ps);
    T = 1.5*(Ps(m).m*Ps(n).T + Ps(n).m*Ps(m).T)/(Ps(m).m + Ps(n).m);
    if isfield(Ps(m), 'v0') && isfield(Ps(n), 'v0')
        T = T + Ps(m).m*Ps(n).m/(Ps(m).m + Ps(n).m)*sum((Ps(m).v0 - Ps(n).v0).^2);
    end
    T = T*e;
    b0 = e^2*abs(Ps(m).q*Ps(n).q)/(4*pi*epsilon_0)./T;
    lnL = log(Ld./b0);
end