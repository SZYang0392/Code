function [Sd, Dd, Pd] = KhiEMCold(P, B, w)
%P = [m; n0; q]
    % epsilon_0 = 8.854187817e-12;
    e = 1.602176565e-19;
    wc = P.q*e.*B./P.m;
    wp2 = 2.899158904791503e-27*P.q.^2.*P.n0./P.m;
    Sd = - wp2./(w.^2 - wc.^2);
    Dd = - Sd.*wc./w;
    Pd = - wp2./w^2;
end