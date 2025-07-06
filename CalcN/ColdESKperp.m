function k_perp = ColdESKperp(w, k_para, B, Ps)
% P : struct with elements q, m, n0.
% w, k_para, B, Ps(*).n0 can be arrays (for safety, of the same size).
    %--------------Determine array size--------------%
    Size = size(w);
    num = numel(w);
    if numel(k_para) > num
        Size = size(k_para);
        num = numel(k_para);
    end
    if numel(B) > num
        Size = size(B);
        num = numel(B);
    end
    if numel(Ps(1).n0) > num
        Size = size(Ps(1).n0);
        num = numel(Ps(1).n0);
    end
    %--------------Calculate S, P--------------%
    Sds = cell(size(Ps));
    Pds = cell(size(Ps));
    S = zeros(Size);
    P = zeros(Size);
    for k = 1 : numel(Ps)
        [Sds{k}, Pds{k}] = Khi(Ps(k), B, w);
        S = S + Sds{k};
        P = P + Pds{k};
    end
    S = S + 1;
    P = P + 1;
    %--------------Calculate k_perp--------------%
    k_perp = sqrt(-P.*k_para.^2./S);
end

function [S1, P1] = Khi(P, B, w)
%P = [m; n0; q]
    % epsilon_0 = 8.854187817e-12;
    e = 1.602176565e-19;
    wc = P.q*e.*B./P.m;
    wp2 = 2.899158904791503e-27*P.q.^2.*P.n0./P.m;
    S1 = - wp2./(w.^2 - wc.^2);
    P1 = - wp2./w.^2;
end