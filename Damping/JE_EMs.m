function JE = JE_EMs(w, k_para, k_perp, B, Ps, E)
    if nargin < 6
        tol_scan = [1e-4, 1e-3, 1e-2, 0.05];
        for tol = tol_scan
            E = Polar(w, k_para, k_perp, B, Ps, 'em', tol);
            if ~isempty(E)
                break;
            end
        end
        
    end
end