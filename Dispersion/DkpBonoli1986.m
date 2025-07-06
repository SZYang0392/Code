function varargout = DkpBonoli1986(w, k_para, k_perp, B, Ps, spec)
% Bonoli P T, Englade R C. 1986, Phys. Fluids, 29: 2937  
% COLLISION EFFECT IS NOT INCLUDED
    % ---------------- Default Input ---------------- %
    if nargin < 6
        Pe = Ps(1);
        Pi = Ps(2:end);
        Alp = false;
    else
        Pe = Ps(1:spec(1));
        Pi = Ps(spec(1)+1:end);
        if numel(spec) >= 2
            Alp = true;
        else
            Alp = false;
        end
    end
    % ---------------- Physical constants ---------------- %
    c = 299792458;
    e = 1.6021766208e-19;
    epsi0 = 8.854187817e-12;
    % ---------------- Parameters ---------------- %
    N_para = k_para.*c./w;
    N_perp = k_perp.*c./w;
    N_para2 = N_para.^2;
    N_perp2 = N_perp.^2;
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
    wpewce = 0;
    for k = 1:numel(Pe)
        wpewce = wpewce + Pe(k).wp2./(Pe(k).wc.^2);
    end
    WPI = 0;
    for k = 1:numel(Pi)
        WPI = WPI + Pi(k).wp2./w2;
    end
    WPE = 0;
    for k = 1:numel(Pe)
        WPE = WPE + Pe(k).wp2./w2;
    end
    % epsilon_perp
    e_perp = 1 + wpewce - WPI;
    % epsilon_para
    e_para = 1 - WPI - WPE;
    % epsilon__xy
    e_xy = 0;
    for k = 1:numel(Pe)
        e_xy = e_xy + Pe(k).wp2./w./Pe(k).wc;
    end
    % ---------------- P0, P2, P4, P6 ---------------- %
    P0 = e_para.*((N_para2 - e_perp).^2 - e_xy.^2);
    P2 = (e_para + e_perp).*(N_para2 - e_perp) + e_xy.^2;
    P4 = e_perp;
    P6E = 0;
    P6I = 0;
    for k = 1:numel(Pe)
        P6E = P6E + Pe(k).wp2 * w2 .* Pe(k).vt2 ./ (Pe(k).wc.^4) /c2;
    end
    for k = 1:numel(Pi)
        P6I = P6I + Pi(k).wp2 .* Pi(k).vt2 ./w2 / c2;
    end
    P6 = - 0.375*P6E - 1.5*P6I;
    Dre = P6.*N_perp2.^3 + P4.*N_perp2.^2 + P2.*N_perp2 + P0;
    % ================ Imaginary Part of the Dispersion relation ================ %
    Dim = 0;
    ktot = sqrt(k_para.^2 + k_perp.^2);
    vphi = w./ktot;
    vphi_para = w./k_para;
    N2 = N_perp2 + N_para.^2;
    DsIm = zeros(size(Ps));
    jIm = 0;
    for k = 1:numel(Pe)
        jIm = jIm + 1;
        DsIm(jIm) = 2*sqrt(pi) * (c2./Pe(k).vt2) .* (Pe(k).wp2./w2) .* N2 .* (vphi_para./sqrt(Pe(k).vt2)) .* exp(-vphi_para.^2./Pe(k).vt2);
        Dim = Dim + DsIm(jIm);
    end
    for k = 1:numel(Pi) - 1
        jIm = jIm + 1;
        DsIm(jIm) = 2*sqrt(pi) * Pi(k).wp2./w2 .* N_perp2.^2 .* (vphi./sqrt(Pi(k).vt2)).^3 .* exp(-vphi.^2./Pi(k).vt2);
        Dim = Dim + DsIm(jIm);
    end
    jIm = jIm + 1;
    k = numel(Pi);
    if Alp
        DsIm(jIm) = N2.^2*ImA(w, Pi(k).n0, Pe(1).T);
        Dim = Dim + DsIm(jIm);
    else
        DsIm(jIm) = 2*sqrt(pi) * Pi(k).wp2./w2 .* N_perp2.^2 .* (vphi./sqrt(Pi(k).vt2)).^3 .* exp(-vphi.^2./Pi(k).vt2);
        Dim = Dim + DsIm(jIm);
    end
    D = Dre + 1i*Dim;
    % ================ Derivative of the Real Part of the Dispersion relation ================ %
    DkpDre = (6*P6.*N_perp2.^2 + 4*P4.*N_perp2 + 2*P2).*N_perp.*(c./w);
    varargout{1} = D;
    varargout{2} = DkpDre;
    varargout{3} = 1i*DsIm;
end

function Di = ImA(w, na, Te_eV)
% Calculate Di caused by alpha particles assuming that nD = nT = 0.5*ne
    e = 1.6021766208e-19;
    epsi0 = 8.854187817e-12;
    ma = 6.642155684e-27;
    v0_vc2 = 3.5e6/32.757585690687492./Te_eV;
    wpa2 = (2*e)^2*na/epsi0./ma;
    Di = -pi*wpa2./(w.^2)./log(0.63*v0_vc2);
end