function Psit = Psip2t(Psi, q)
% Calculate toroidal flux with poloidal flux and safety factor.
    Psit = nan(size(Psi));
    Psit(1) = 0;
    qmid = 0.5*(q(1:end-1) + q(2:end));
    dPsip = 2*pi*(Psi(2:end) - Psi(1:end-1));
    dPsit = qmid.*dPsip;
    for j = 2:numel(Psit)
        Psit(j) = Psit(j-1) + dPsit(j-1);
    end
end