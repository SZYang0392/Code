function Khi = KhiMaxwell_sum(w, k_para, k_perp, B, Ps, Psi)
% Calculate susceptibility tensor for multiple P species.
    if nargin < 6
        Psi = 0;
    end
    Khi = zeros(3, 3) + 1i*zeros(3, 3);
    for k = 1:numel(Ps)
        Khi = Khi + KhiMaxwell(w, k_para, k_perp, B, Ps(k), Psi);
    end
end