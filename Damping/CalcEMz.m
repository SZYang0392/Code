function EMz = CalcEMz(w, k_para, B, Ps, Opts)
%% This function calculate the following properties with EM fzero method in order:
%   Re(k_perp)
%   Khih(wr, kr), d[Khih(wr, kr)]/dw, d[Khih(wr, kr)]/dk_perp, Khia(wr, kr)
%   Khi, Khi_s
%   Epsilonh(wr, kr), d[Epsilonh(wr, kr)]/dw, d[Epsilonh(wr, kr)]/dk_perp, Epsilona(wr, kr)
%   Epsilon, d[wr*Epsilonh(wr, kr)]/dw
%   Re[D(wr, kr)], dRe[D(wr, kr)]/dw, dRe[D(wr, kr)]/dk_perp, Im[D(wr, kr)]
%   D(w, k) = Re[D(wr, kr)] + iwi*dRe[D(wr, kr)]/dw + iki*dRe[D(wr, kr)]/dk_perp + Im[D(wr, kr)]
%   wi, ki
%   vg = -  dRe(D(wr, kr))/dk_perp  /  dRe(D(wr, kr))/dw
%   E (D(w, k)*E = 0), <J.E>

%% Check arguments
arguments
    w (1, 1) {mustBeReal, mustBePositive}
    k_para (1, 1) {mustBeReal}
    B  (1, 1) {mustBeReal, mustBePositive}
    Ps 
    Opts.Nscan (1, 1) {mustBePositive, mustBeInteger} = 11
    Opts.range (1, 2) {mustBeReal, mustBeNonNan} = [0.8, 3.5]
    Opts.dw (1, 1) {mustBeReal, mustBePositive} = 500
    Opts.dw2w (1, 1) {mustBeReal, mustBePositive} = 0.001
    Opts.dk (1, 1) {mustBeReal, mustBePositive} = 10
    Opts.dk2k (1, 1) {mustBeReal, mustBePositive} = 0.001
    Opts.k_perp_c
end
    
    %% Initialize dw, k_perp_c and P
    dw = min(Opts.dw, Opts.dw2w*w);
    if isfield(Opts, 'k_perp_EM')
        k_perp_c = Opts.k_perp_c;
    else
        k_perp_c = ColdEMKperp(w, k_para, B, Ps);
    end
    P = Ps;
    for k = 1:numel(P)
        P(k).T = P(k).T(1);
        P(k).n0 = P(k).n0(1);
    end
    nr = numel(Ps(1).n0);
    np = numel(Ps);
    P = repmat(P, nr, 1);
    for j = 1:nr
        for k = 1:numel(Ps)
            P(j, k).T = Ps(k).T(j);
            P(j, k).n0 = Ps(k).n0(j);
        end
    end
    %% Initialize output
    EMz.Cutline0 = '---------------------------- Parameters ----------------------------';
    EMz.w = w;
    EMz.k_para = k_para;
    EMz.B = B;
    EMz.Ps = Ps;
    EMz.dw = dw;
    EMz.Cutline1 = '---------------------------- Propagating and Damping ----------------------------';
    EMz.k_perp = nan(2, nr);
    EMz.ki_perp = nan(2, nr);
    EMz.gamma = nan(2, nr);
    EMz.vg = nan(2, nr);
    EMz.Cutline2 = '---------------------------- Dielectric Properties ----------------------------';
    nanC33 = repmat(nan + nan*1i, 3, 3, nr);
    zero33 = repmat(0 + 0*1i, 3, 3, nr);
    eye33 = repmat(eye(3), 1, 1, nr);
    % D(w, k)
    EMz.ReD = nan(2, nr);
    EMz.ReDw = nan(2, nr);
    EMz.ReDk = nan(2, nr);
    EMz.ImD = nan(2, nr);
    EMz.D = nan(2, nr);
    % Khi_s
    EMz.Khih_s = cell(2, np);
    EMz.Khihw_s = cell(2, np);
    EMz.Khihk_s = cell(2, np);
    EMz.Khia_s = cell(2, np);
    EMz.Khi_s = cell(2, np);
    for j = 1:2
        for l = 1:np
            EMz.Khih_s{j, l} = nanC33;
            EMz.Khihw_s{j, l} = nanC33;
            EMz.Khihk_s{j, l} = nanC33;
            EMz.Khia_s{j, l} = nanC33;
            EMz.Khi_s{j, l} = nanC33;
        end
    end
    % Khi
    EMz.Khih = {zero33; zero33};
    EMz.Khihw = {zero33; zero33};
    EMz.Khihk = {zero33; zero33};
    EMz.Khia = {zero33; zero33};
    EMz.Khi = {zero33; zero33};
    % Epsilon
    EMz.Epsh = {eye33; eye33};
    EMz.Epshw = {zero33; zero33};
    EMz.Epshk = {zero33; zero33};
    EMz.Epsa = {zero33; zero33};
    EMz.Eps = {eye33; eye33};
    EMz.dwEpshdw = {zero33; zero33};
    % D3x3
    EMz.D33 = EMz.Eps;
    EMz.Cutline3 = '---------------------------- Electric field ----------------------------';
    EMz.E = repmat(nan + 1i*nan, 3, nr, 2);
    EMz.JEs = nan(np, nr, 2);
    EMz.JE = zeros(2, nr);
    %% Calculate Re(k_perp) and gather
    k_perp_EM0 = nan(Opts.Nscan, nr);
    for j = 1:nr
        % Set range of scan
        if real(k_perp_c(2, j)) > 0 && real(k_perp_c(1, j)) > 0
            kscan = linspace(Opts.range(1)*real(k_perp_c(2, j)), Opts.range(2)*real(k_perp_c(1, j)), Opts.Nscan);
        elseif real(k_perp_c(2, j)) <= 0 && real(k_perp_c(1, j)) > 0
            kscan = linspace(Opts.range(1)*real(k_perp_c(1, j)), Opts.range(2)*real(k_perp_c(1, j)), Opts.Nscan);
        else
            continue;
        end
        % Scan k_perp to find the solution
        Dscan = nan(size(kscan));
        f = @(k_perpt) real(DEMMaxwell(w, k_para, k_perpt, B, P(j, :)));
        Dscan(1) = f(kscan(1));
        for l = 2:Opts.Nscan
            Dscan(l) = f(kscan(l));
            if Dscan(l-1)*Dscan(l) <= 0
                ki_perp_EMs1 = fzero(f, [kscan(l-1), kscan(l)]);
                k_perp_EM0(l, j) = ki_perp_EMs1;
            end
        end
    end
    EMz.k_perp = Gather(k_perp_EM0);
    %% Calculate damping
    for l = 2:-1:1
        for j = 1:nr
            k_perp = EMz.k_perp(l, j);
            dk = min(Opts.dk, Opts.dk2k*k_perp);
            D = DEMMaxwell(w, k_para, k_perp, B, P(j, :));
            EMz.ReD(l, j) = real(D);
            EMz.ImD(l, j) = imag(D);
            EMz.ReDw(l, j) = real(DEMMaxwell(w+dw, k_para, k_perp, B, P(j, :)) - DEMMaxwell(w-dw, k_para, k_perp, B, P(j, :)))/2/dw;
            EMz.ReDk(l, j) = real(DEMMaxwell(w, k_para, k_perp+dk, B, P(j, :)) - DEMMaxwell(w, k_para, k_perp-dk, B, P(j, :)))/2/dk;
        end
    end
    EMz.gamma = -EMz.ImD./EMz.ReDw;
    EMz.ki_perp = -EMz.ImD./EMz.ReDk;
    EMz.vg = -EMz.ReDk./EMz.ReDw;
    EMz.D = EMz.ReD(l, j) + 1i.*EMz.gamma.*EMz.ReDw + 1i.*EMz.ki_perp.*EMz.ReDk + 1i*EMz.ImD;
    %% Calculate Khi_s
    for l = 2:-1:1
        for j = 1:nr
            if isnan(EMz.k_perp(l, j))
                continue;
            end
            k_perp = EMz.k_perp(l, j);
            dk = min(Opts.dk, Opts.dk2k*k_perp);
            for k = 1:np
                EMz.Khih_s{l, k}(:, :, j) = KhiMaxwell(w, k_para, k_perp, B, P(j, k));
                EMz.Khia_s{l, k}(:, :, j) = 0.5*(EMz.Khih_s{l, k}(:, :, j) - EMz.Khih_s{l, k}(:, :, j)');
                EMz.Khih_s{l, k}(:, :, j) = 0.5*(EMz.Khih_s{l, k}(:, :, j) + EMz.Khih_s{l, k}(:, :, j)');
                EMz.Khihw_s{l, k}(:, :, j) = (KhiMaxwell(w+dw, k_para, k_perp, B, P(j, k)) - KhiMaxwell(w-dw, k_para, k_perp, B, P(j, k)))/2/dw;
                EMz.Khihw_s{l, k}(:, :, j) = 0.5*(EMz.Khihw_s{l, k}(:, :, j) + EMz.Khihw_s{l, k}(:, :, j)');
                EMz.Khihk_s{l, k}(:, :, j) = (KhiMaxwell(w, k_para, k_perp+dk, B, P(j, k)) - KhiMaxwell(w, k_para, k_perp-dk, B, P(j, k)))/2/dk;
                EMz.Khihk_s{l, k}(:, :, j) = 0.5*(EMz.Khihk_s{l, k}(:, :, j) + EMz.Khihk_s{l, k}(:, :, j)');
                EMz.Khi_s{l, k}(:, :, j) = EMz.Khih_s{l, k}(:, :, j) + 1i*EMz.gamma(l, j)*EMz.Khihw_s{l, k}(:, :, j) + 1i*EMz.ki_perp(l, j)*EMz.Khihk_s{l, k}(:, :, j) + EMz.Khia_s{l, k}(:, :, j);
            end
        end
    end
    %% Calculate Khi and Epsilon
    gamma3 = reshape(EMz.gamma, 2, 1, []);
    ki3 = reshape(EMz.ki_perp, 2, 1, []);
    for l = 2:-1:1
        for k = 1:np
            EMz.Khih{l} = EMz.Khih{l} + EMz.Khih_s{l, k};
            EMz.Khihw{l} = EMz.Khihw{l} + EMz.Khihw_s{l, k};
            EMz.Khihk{l} = EMz.Khihk{l} + EMz.Khihk_s{l, k};
            EMz.Khia{l} = EMz.Khia{l} + EMz.Khia_s{l, k};
        end
        EMz.Khi{l} = EMz.Khih{l} + 1i.*gamma3(l, :, :).*EMz.Khihw{l} + 1i.*ki3(l, :, :).* EMz.Khihk{l} + EMz.Khia{l};
        % Epsilon
        EMz.Epsh{l} = EMz.Epsh{l} + EMz.Khih{l};
        EMz.Epsa{l} = EMz.Epsa{l} + EMz.Khia{l};
        EMz.Eps{l} = EMz.Eps{l} + EMz.Khi{l};
        EMz.dwEpshdw{l} = EMz.dwEpshdw{l} + EMz.Epsh{l} + w*EMz.Khihw{l};
    end
    EMz.Epshw = EMz.Khihw;
    EMz.Epshk = EMz.Khihk;
    c = 299792458;
    EMz.D33 = EMz.Eps;
    for l = 2:-1:1
        for j = 1:nr
            k_perp = EMz.k_perp(l, j);
            N = [k_perp + 1i*EMz.ki_perp(l, j), 0, k_para]*c/w;
            EMz.D33{l}(:, :, j) = N.'*N - N*N.'*eye(3) + EMz.Khi{l}(:, :, j);
        end
    end
    %% Calculate E and J.E
    epsilon_0 = 8.854187817e-12;
    for l = 2:-1:1
        for j = 1:nr
            E0 = null(EMz.D33{l}(:, :, j));
            if isempty(E0)
                % EMz.E = repmat(nan + 1i*nan, 3, nr, 2);
                % EMz.JEs = nan(np, nr, 2);
                % EMz.JE = nan(2, nr);
                EMz.E(:, j, l) = nan;
                EMz.JEs(:, j, l) = nan;
                EMz.JE(l, j) = nan;
            else
                E0 = E0(:, 1);
                EMz.E(:, j, l) = E0;
                for k = 1:np
                    sigmah = -1i*epsilon_0*w*EMz.Khia_s{l, k}(:, :, j);
                    EMz.JEs(k, j, l) = real(E0'*sigmah*E0);
                    EMz.JE(l, j) = EMz.JE(l, j) + EMz.JEs(k, j, l);
                end
            end
        end
    end
end

%% Functions
% Function to gather slow and fast wave solutions.
function ki_G = Gather(ki_n)
    if isreal(ki_n(1))
        ki_G = nan(2, size(ki_n, 2), size(ki_n, 3));
    else
        ki_G = repmat(nan + 1i*nan, [2, size(ki_n, 2), size(ki_n, 3)]);
    end
    for k = 1:size(ki_n, 3)
        Check = ~isnan(ki_n(:, :, k));
        Num = sum(Check, 1);
        Nummax = max(Num, [], 'all');
        ki_n2_0 = nan(max(Nummax, 2), size(ki_n, 2));
        for l = 1:Nummax
            Ind = Num == l;
            if ~isempty(Ind)
                ki_n2_1 = ki_n(:, Ind, k);
                Ind1 = ~isnan(ki_n2_1);
                ki_n2_2 = reshape(ki_n2_1(Ind1), l, sum(Ind));
                if l == 1
                    ki_n2_0(2, Ind) = ki_n2_2;
                else
                    ki_n2_0(1:l, Ind) = ki_n2_2;
                end
            end
        end
        ki_G(:, :, k) = ki_n2_0(1:2, :);
    end
end