function varargout = DEMCold(w, k_para, k_perp, B, Ps)
% P : struct with elements q, m, n0.
% w, k_para, B, Ps(*).n0 can be arrays (for safety, of the same size).
    % addpath('../Dispersion');
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
    %--------------Refractive Index--------------%
    c = 299792458;
    N_para = k_para.*c./w;
    N_perp = k_perp.*c./w;
    %--------------Calculate S, D, P--------------%
    Sds = cell(size(Ps));
    Dds = cell(size(Ps));
    Pds = cell(size(Ps));
    S = zeros(Size);
    D = zeros(Size);
    P = zeros(Size);
    for k = 1 : numel(Ps)
        [Sds{k}, Dds{k}, Pds{k}] = KhiEMCold(Ps(k), B, w);
        S = S + Sds{k};
        D = D + Dds{k};
        P = P + Pds{k};
    end
    S = S + 1;
    P = P + 1;
    %--------------Calculate P4, P2, P0--------------%
    P4 = S;
    P2 = (N_para.^2 - S).*(S + P) + D.^2;
    P0 = P.*((N_para.^2 - S).^2 - D.^2);

    %====================Calculate D====================%
    Dispersion = P4.*N_perp.^4 + P2.*N_perp.^2 + P0;
    varargout{1} = Dispersion;
    varargout{2} = [P0, P2, P4];
    varargout{3} = [S, D, P];
end