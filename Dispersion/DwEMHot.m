function varargout = DwEMHot(w, k_para, k_perp, B, Ps)
% Calculate EM Dispersion relation with Maxwellian distribution.
% w, N_para, N_perp, B, Ps(k).* are scalars.
    c = 299792458;
    N_para = k_para.*c./w;
    N_perp = k_perp.*c./w;
    N = [N_perp     0     N_para];
    Khi = zeros(3, 3, numel(Ps));
    Khiw = zeros(3, 3, numel(Ps));
    for k = 1:numel(Ps)
        [Khi(:, :, k), Khiw(:, :, k)] = KhiwMaxwell(w, k_para, k_perp, B, Ps(k), 0);
    end
    % D
    Khi_Sum = sum(Khi, 3);
    Dielectric = eye(3) + Khi_Sum;
    DM = N.'*N - eye(3).*(N_para.^2 + N_perp.^2) + Dielectric;
    D = det(DM);
    % dD/dw
    Khiw_Sum = sum(Khiw, 3);
    DwDM = (N.'*N - eye(3).*(N_para.^2 + N_perp.^2)).*(-2./w) + Khiw_Sum;
    DwDM1 = [DwDM(1, :); DM([2, 3], :)];
    DwDM2 = [DM(1, :); DwDM(2, :); DM(3, :)];
    DwDM3 = [DM([1, 2], :); DwDM(3, :)];
    DwD = det(DwDM1) + det(DwDM2) + det(DwDM3);
    varargout{1} = D;
    varargout{2} = DwD;
end