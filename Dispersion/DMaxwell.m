function Dielectric = DMaxwell(w, k_para, k_perp, B, Ps, Psi)
% Calculate dielectric tensor.
    if nargin < 6
        Psi = 0;
    end
    Dielectric = eye(3) + KhiMaxwell_sum(w, k_para, k_perp, B, Ps, Psi);
end