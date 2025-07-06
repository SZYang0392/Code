function varargout = DEMMaxwell(w, k_para, k_perp, B, Ps)
% Calculate EM Dispersion relation with Maxwellian distribution.
% w, N_para, N_perp, B, Ps(k).* are scalars.
    c = 299792458;
    N_para = k_para.*c./w;
    N_perp = k_perp.*c./w;
    N = [N_perp, 0, N_para];
    Khi = zeros(3, 3, numel(Ps));
    for k = 1:numel(Ps)
        if isfield(Ps(k), 'mag')
            if Ps(k).mag == 0
                Khi(:, :, k) = KhiESMaxwell(w, k_para, k_perp, Ps(k));
            else
                Khi(:, :, k) = KhiMaxwell(w, k_para, k_perp, B, Ps(k), 0);
            end
        else
            Khi(:, :, k) = KhiMaxwell(w, k_para, k_perp, B, Ps(k), 0);
        end
    end
    Khi_Sum = sum(Khi, 3);
    Dielectric = eye(3) + Khi_Sum;
    Dtensor = N.'*N - eye(3).*(N_para.^2 + N_perp.^2) + Dielectric;
    D = det(Dtensor);

    varargout{1} = D;
    varargout{2} = Khi_Sum;
    varargout{3} = Khi;
    varargout{4} = Dielectric;
    varargout{5} = Dtensor;
end