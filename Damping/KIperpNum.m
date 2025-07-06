function varargout = KIperpNum(w, k_para, k_perp, B, Ps, f, ifspec)
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
    dk = k_perp*1e-4;
    DkD = (f(w, k_para, k_perp + dk, B, Ps) - f(w, k_para, k_perp - dk, B, Ps));
    DkD = DkD/2.0/dk;

    %Damping rate
    KI = imag(D)/real(DkD);
    KIs = imag(Ds)/real(DkD);

    %Output
    varargout{1} = KI;
    varargout{2} = D;
    varargout{3} = KIs;
    varargout{4} = Ds;
end