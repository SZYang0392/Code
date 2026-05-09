function D = DESCold(w, k_para, k_perp, B, Ps)
% P : struct with elements q, m, n0.
% w, k_para, B, Ps(*).n0 can be arrays (for safety, of the same size).
    % addpath('../Dispersion');
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
    S = zeros(Size);
    P = zeros(Size);
    for k = 1 : numel(Ps)
        [Sds{k}, Pds{k}] = KhiESCold(Ps(k), B, w);
        S = S + Sds{k};
        P = P + Pds{k};
    end
    S = S + 1;
    P = P + 1;
    %--------------Calculate Dispersion--------------%
    D = (c./w).^2.*(S.*k_perp.^2 + P.*k_para.^2);
end