function lnL = LnLambda(Ps, m, n)
% Calculate Coulomb logarithm
    if nargin <= 2
        P1 = Ps;
        P2 = m;
    else
        P1 = Ps(m);
        P2 = Ps(n);
    end
    e = 1.602176565e-19;
    epsilon_0 = 8.854187817e-12;

    Ld = Debye(Ps);
    T = 1.5*(P1.m*P2.T + P2.m*P1.T)/(P1.m + P2.m);
    if isfield(P1, 'v0') && isfield(P2, 'v0')
        T = T + P1.m*P2.m/(P1.m + P2.m)*sum((P1.v0 - P2.v0).^2);
    end
    T = T*e;
    b0 = e^2*abs(P1.q*P2.q)/(4*pi*epsilon_0)./T;
    lnL = log(Ld./b0);
end