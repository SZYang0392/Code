function varargout = DwESHot(w, k_para, k_perp, B, Ps)
% Calculate ES Dispersion relation with Maxwellian distribution.
% w, N_para, N_perp, B, Ps(k).* are scalars.
    c = 299792458;
    N_para = k_para.*c./w;
    N_perp = k_perp.*c./w;
    N = [N_perp, 0, N_para];
    Khi = zeros(3, 3, numel(Ps));
    Khiw = zeros(3, 3, numel(Ps));
    Ds = zeros(1, numel(Ps));
    for k = 1:numel(Ps)
        [Khi(:, :, k), Khiw(:, :, k)] = KhiwMaxwell(w, k_para, k_perp, B, Ps(k), 0);
        Ds(k) = N*Khi(:, :, k)*(N.');
    end
    % D
    D = sum(Ds, 'all') + N*eye(3)*(N.');
    varargout{1} = D;
    % dD/dW
    Khiw_Sum = sum(Khiw, 3);
    DwD = N*Khiw_Sum*(N.') - (2./w).*D;
    varargout{2} = DwD;
    % Species
    varargout{3} = Ds.';
end