function varargout = DLandau(w, k, Ps)
    % ES dispersion for nonmagnetized plasma.
    e = 1.602176565e-19;
    epsilon_0 = 8.854187817e-12;

    vphi = w/k;
    D = 0;
    for j = 1:numel(Ps)
        if ~isfield(Ps(j), 'vt')
            Ps(j).vt = sqrt(2*Ps(j).T*e/Ps(j).m);
        elseif isempty(Ps(j).vt)
            Ps(j).vt = sqrt(2*Ps(j).T*e/Ps(j).m);
        end
        if ~isfield(Ps(j), 'wp')
            Ps(j).wp = abs(Ps(j).q)*e*sqrt(Ps(j).n0/Ps(j).m/epsilon_0);
        elseif isempty(Ps(j).wp)
            Ps(j).wp = abs(Ps(j).q)*e*sqrt(Ps(j).n0/Ps(j).m/epsilon_0);
        end

        Ps(j).zeta = vphi/Ps(j).vt;
        Ps(j).Khi = 2*(Ps(j).wp/k/Ps(j).vt)^2*(1 + Ps(j).zeta*Z_plasma(Ps(j).zeta));
        D = D + Ps(j).Khi;
    end

    D = D + 1;

    varargout{1} = D;
    varargout{2} = Ps;
end