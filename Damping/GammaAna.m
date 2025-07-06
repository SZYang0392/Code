function varargout = GammaNum(w, k_para, k_perp, B, Ps, f, ifspec)
%Calculate damping rate gamma numerically.
%f : function handle of dispersion relation.
    if nargin < 7
        ifspec = false;
    end

    %Calculate Dispersion and Direvative
    if ifspec
        [D, DwD, Ds] = f(w, k_para, k_perp, B, Ps);
    else
        [D, DwD] = f(w, k_para, k_perp, B, Ps);
        Ds = [];
    end

    %Damping rate
    Gamma = imag(D)/real(DwD);
    Gammas = imag(Ds)/real(DwD);

    %Output
    varargout{1} = Gamma;
    varargout{2} = D;
    varargout{3} = Gammas;
    varargout{4} = Ds;
end