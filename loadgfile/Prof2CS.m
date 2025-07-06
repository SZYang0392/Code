function f_CS = Prof2CS(Psi_prof, f_prof, Psi_cs, inLCFS)
% Calculate field on the cross section with values along the minor radius.
    [nz, nr] = size(Psi_cs);
    Psi_cs_1D = reshape(Psi_cs, 1, []);
    f_CS_1D = spline(Psi_prof, f_prof, Psi_cs_1D);
    f_CS = reshape(f_CS_1D, nz, nr);
    if nargin > 3
        f_CS(~inLCFS) = nan;
    end
end