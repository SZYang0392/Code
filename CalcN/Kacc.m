function k_acc = Kacc(w, B, Ps)
% ZHAI Xuemei, 2019, Ph.D. Thesis
% Calculate accessable condition of LHW
% N_para^2     >     S + 2D*sqrt(S / -P) - D^2/P     =     (sqrt(S) + D / sqrt(-P))^2
% P : struct with fields q, m, n0.
% w, k_para, B, Ps(*).n0 can be arrays of the same size.
    %--------------Refractive Index--------------%
    c = 299792458;
    %--------------Calculate S, D, P--------------%
    Sds = cell(size(Ps));
    Dds = cell(size(Ps));
    Pds = cell(size(Ps));
    S = 1;
    D = 0;
    P = 1;
    for k = 1 : numel(Ps)
        [Sds{k}, Dds{k}, Pds{k}] = KhiEMCold(Ps(k), B, w);
        S = S + Sds{k};
        D = D + Dds{k};
        P = P + Pds{k};
    end
    %--------------Calculate k_para at cut off--------------%
    N_acc = abs(sqrt(S) + D./sqrt(-P));
    k_acc = N_acc.*w./c;
end

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