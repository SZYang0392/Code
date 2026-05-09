function [S1, P1] = KhiESCold(P, B, w)
%P = [m; n0; q]
    % epsilon_0 = 8.854187817e-12;
    e = 1.602176565e-19;
    wc = P.q*e.*B./P.m;
    wp2 = 2.899158904791503e-27*P.q.^2.*P.n0./P.m;
    S1 = - wp2./(w.^2 - wc.^2);
    P1 = - wp2./w.^2;
end