function varargout = DwESCold(w, k_para, k_perp, B, Ps)
% P : struct with elements q, m, n0.
% w, k_para, B, Ps(*).n0 can be arrays (for safety, of the same size).
    c = 299792458;
    %--------------Determine array size--------------%
    Size = size(w);
    num = numel(w);
    if numel(k_para) > num
        Size = size(k_para);
        num = numel(k_para);
    end
    if numel(k_perp) > num
        Size = size(k_perp);
        num = numel(k_perp);
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
    DwSds = cell(size(Ps));
    DwPds = cell(size(Ps));
    Ds = zeros(numel(Ps), num);
    S = zeros(Size);
    P = zeros(Size);
    DwD = zeros(Size);
    for k = 1 : numel(Ps)
        [Sds{k}, Pds{k}, DwSds{k}, DwPds{k}] = Khiw(Ps(k), B, w);
        S = S + Sds{k};
        P = P + Pds{k};
        DwD = DwD + DwSds{k}.*k_perp.^2 + DwPds{k}.^k_para.^2;
        Ds(k, :) = Sds{k}.*k_perp.^2 + Pds{k}.^k_para.^2;
    end
    S = S + 1;
    P = P + 1;
    %--------------Calculate Dispersion--------------%
    D = (c./w).^2.*(S.*k_perp.^2 + P.*k_para.^2);
    Ds = (c./w).^2.*Ds;
    DwD = DwD.*(c./w).^2 - (2./w).*D;

    varargout{1} = D;
    varargout{2} = DwD;
    varargout{3} = Ds;
end

function [S1, P1, DwS1, DwP1] = Khiw(P, B, w)
%P = [m; n0; q]
    % epsilon_0 = 8.854187817e-12;
    e = 1.602176565e-19;
    wc = P.q*e.*B./P.m;
    wp2 = 2.899158904791503e-27*P.q.^2.*P.n0./P.m;
    S1 = - wp2./(w.^2 - wc.^2);
    P1 = - wp2./w.^2;
    DwS1 = wp2./(w.^2 - wc.^2).^2.*2.*w;
    DwP1 = -(2./w).*P1;
end