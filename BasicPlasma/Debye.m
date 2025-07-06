function varargout = Debye(Ps)
% Calculate Debye Length
    e = 1.602176565e-19;
    epsilon_0 = 8.854187817e-12;

    L = zeros(1, numel(Ps(1).n0));
    Ls = nan(numel(Ps), numel(Ps(1).n0));
    for k = 1:numel(Ps)
        Ls(k, :) = sqrt(epsilon_0*Ps(k).T./(e*Ps(k).q^2.*Ps(k).n0));
        L = L + 1./(Ls(k, :).^2);
    end
    L = 1./sqrt(L);

    if nargout == 1
        varargout{1} = Ls(1, :);
    else
        varargout{1} = L;
        varargout{2} = Ls;
    end
end