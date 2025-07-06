function varargout = GammaNum(w, k_para, k_perp, B, Ps, f, ifspec)
%Calculate damping rate gamma numerically.
%f : function handle of dispersion relation.
    if nargin < 7
        ifspec = false;
    end

    %Calculate Dispersion
    if ifspec
        [D, Ds] = f(w, k_para, k_perp, B, Ps);
    else
        D = f(w, k_para, k_perp, B, Ps);
        Ds = [];
    end

    %Calculate dirivative
    dw = 300;
    DwD = (f(w + dw, k_para, k_perp, B, Ps) - f(w - dw, k_para, k_perp, B, Ps));
    DwD = DwD/2.0/dw;

    %Damping rate
    Gamma = imag(D)/real(DwD);
    Gammas = imag(Ds)/real(DwD);

    %Output
    varargout{1} = Gamma;
    varargout{2} = D;
    varargout{3} = Gammas;
    varargout{4} = Ds;
end