function k_perp = ColdEMKperp(w, k_para, B, Ps)
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

    %====================Calculate N_perp====================%
    %--------------------Analytically--------------------%
    N1 = - P2./2./P4;
    N2 = sqrt(P2.^2 - 4*P0.*P4)./2./P4;
    N_fast = sqrt(N1 - N2);
    N_slow = sqrt(N1 + N2);
    %--------------------Numerically--------------------%
    % N_perp = roots([P4, P2, P0]);
    % N_perp = sqrt(sort(N_perp, 'ComparisonMethod', 'real'));

    %====================Put slow wave at N_perp(1, :)====================%
    % N_perp.^2 of the fast wave is smaller than that of the slow wave.
    % The first element of N_perp denotes the perpendicular refraction index of the slow wave.
    ind = abs(real(N_slow)) < abs(real(N_fast));
    N3 = zeros(size(N1));
    N3(ind) = N_slow(ind);
    N_slow(ind) = N_fast(ind);
    N_fast(ind) = N3(ind);
    N_perp = [N_slow; N_fast];

    %--------------Transfer N_perp to k_perp--------------%
    k_perp = N_perp.*w./c;
end