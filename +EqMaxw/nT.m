function varargout = nT(Ps, Ph, E0)
    % E0 : maximum energy in eV
    if nargin < 3
        if isfield(Ph, 'u0')
            if ~isempty(Ph.u0)
                E0 = 0.5*Ph.m*Ph.u0^2;
            else
                fprintf('In EqMaxw/Temp.m : Ph.u0 is empty\n');
            end
        else
            fprintf('In EqMaxw/Temp.m : Error inputting E0\n');
        end
    end
    % Ec : critical energy in eV, PATH 'SlowingDown' needed.
    Ec0 = Ec(Ph.m, Ps);
    Tsd1 = EqMaxw.Tsd(E0, Ec0);
    Etasd1 = Ph.n0./Ps(1).n0;
    Etasd2 = EqMaxw.Eta_q(Etasd1, E0, Tsd1);
    T = EqMaxw.Tsd_q(E0, Tsd1, Etasd1, Etasd2);
    n = Ps(1).n0.*Etasd2;
    varargout{1} = n;
    varargout{2} = T;
    varargout{3} = Etasd1;
    varargout{4} = Etasd2;
    varargout{5} = Tsd1;
end