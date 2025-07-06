function E = Eparticle(fnum, tscale)
    % Load (particle energy / m) from "particle*****.dat"
    % fnum : specify files to be loaded
    % tscale : fnum*tscale = wci*t; tscale = dt*me/mi*ndiagp
    
    if nargin < 2
        tscale = 1;
    end

    E.stime = fnum*tscale;
    P = loadp_f(fnum(1));
    nele = numel(P.ele);
    nion = numel(P.ion);
    E.Eele = zeros(nele, numel(fnum));
    E.Eion = zeros(nion, numel(fnum));

    for k = 1:numel(fnum)
        P = loadp_f(fnum(k));
        for j = 1:nele
            E.Eele(j, k) = 1.5*mean(P.ele(j).v.^2, 'all');
        end
        for j = 1:nion
            E.Eion(j, k) = 1.5*mean(P.ion(j).v.^2, 'all');
        end
        fprintf('%d / %d\n', k, numel(fnum));
    end
end