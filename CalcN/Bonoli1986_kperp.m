function k_perp = Bonoli1986_kperp(w, k_para, B, Pe, Pi)
% Bonoli P T, Englade R C. 1986, Phys. Fluids, 29: 2937  
    % ---------------- Default Input ---------------- %
    if nargin < 5
        Pi = Pe(2:end);
        Pe = Pe(1);
    end
    % ---------------- Physical constants ---------------- %
    c = 299792458;
    e = 1.6021766208e-19;
    epsi0 = 8.854187817e-12;
    % ---------------- Parameters ---------------- %
    N_para = k_para.*c./w;
    w2 = w.^2;
    c2 = c^2;
    for k = 1:numel(Pe)
        Pe(k).wc = Pe(k).q.*e.*B./Pe(k).m;
        Pe(k).wp2 = (Pe(k).q.*e).^2.*Pe(k).n0./epsi0./Pe(k).m;
        Pe(k).vt2 = 2*Pe(k).T*e./Pe(k).m;
    end
    for k = 1:numel(Pi)
        Pi(k).wc = Pi(k).q.*e.*B./Pi(k).m;
        Pi(k).wp2 = (Pi(k).q.*e).^2.*Pi(k).n0./epsi0./Pi(k).m;
        Pi(k).vt2 = 2*Pi(k).T*e./Pi(k).m;
    end
    % ================ Real Part of the Dispersion relation ================ %
    % epsilon_perp
    e_perp = 1;
    for k = 1:numel(Pe)
        e_perp = e_perp + Pe(k).wp2./(Pe(k).wc.^2);
    end
    for k = 1:numel(Pi)
        e_perp = e_perp - Pi(k).wp2./w2;
    end
    % epsilon_para
    e_para = 1;
    for k = 1:numel(Pe)
        e_para = e_para - Pe(k).wp2./w2;
    end
    for k = 1:numel(Pi)
        e_para = e_para - Pi(k).wp2./w2;
    end
    % epsilon__xy
    e_xy = 0;
    for k = 1:numel(Pe)
        e_xy = e_xy + Pe(k).wp2./w./Pe(k).wc;
    end
    % ---------------- P0, P2, P4, P6 ---------------- %
    P0 = e_para.*((N_para.^2 - e_perp).^2 - e_xy.^2);
    P2 = (e_para + e_perp).*(N_para.^2 - e_perp) + e_xy.^2;
    P4 = e_perp;
    P6 = 0;
    for k = 1:numel(Pe)
        P6 = P6 - 0.375 * Pe(k).wp2 * w2 .* Pe(k).vt2 ./ (Pe(k).wc.^4) /c2;
    end
    for k = 1:numel(Pi)
        P6 = P6 - 1.5 * Pi(k).wp2 .* Pi(k).vt2 ./w2 / c2;
    end
    % ================ Solution ================ %
    N_perp2 = roots([P6, P4, P2, P0]);
    k_perp = sqrt(N_perp2)*w/c;
end