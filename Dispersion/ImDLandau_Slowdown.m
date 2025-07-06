function ImD = ImDLandau_Slowdown(w, k, P)
    % Calculate ImD for slowing-down distribution for Landau damping.
    % -------------------- Physical Constants -------------------- %
    c = 299792458;
    e = 1.602176565e-19;
    epsilon_0 = 8.854187817e-12;

    % -------------------- Some Factors -------------------- %
    ktot2 = k.^2;
    N2 = ktot2.*(c/w).^2;
    Fac = 3./4./log(1 + (P.u0./P.uc).^3);
    wp2 = (P.q.*e).^2.*P.n0./epsilon_0./P.m;
    ImD = N2.*Fac.*(wp2./k.^2).*(w./k)./(P.uc.^3 + (w./k).^3);
end