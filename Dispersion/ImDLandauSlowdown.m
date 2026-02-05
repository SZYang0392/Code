function ImD = ImDLandauSlowdown(w, ktot, P, k_para, ifautomatic)
% Calculate unmagnetized Landau damping with slowing-down distribution
% ImD = N^2 N.\Chi^a.N
    errmax = 5e-3;
    if nargin > 3
        k_perp = ktot;
        ktot = sqrt(k_para.^2 + k_perp.^2);
    end
    if nargin < 5
        ifautomatic = false;
    end
    % -------------------- Physical Constants -------------------- %
    c = 299792458;
    e = 1.602176565e-19;
    epsilon_0 = 8.854187817e-12;

    % -------------------- Some Key Factors -------------------- %
    N = ktot*c./w;
    wp2 = (P.q*e).^2.*P.n0./epsilon_0./P.m;
    w2 = w.^2;
    uc = P.uc./P.u0;
    zeta = w./ktot./P.u0;

    isotropic = ~isfield(P, 'DLambda') || P.DLambda == 0 || isnan(P.DLambda) || isempty(P.DLambda);
    if isotropic
        % Vector calculation is available in isotropic case
        ImD = N.^4.*wp2./w2.*(1.5*pi)./log(1 + 1./(uc.^3)).*zeta.^3.*(1./(zeta.^3 + uc.^3) - 1./(1 + uc.^3));
        Ind = zeta > 1;
        ImD(Ind) = 0;
    else
        % Calculate B
        if ifautomatic
            ftheta = @(x) exp(-((sin(x).^2 - P.Lambda0)./P.DLambda).^2).*sin(x);
            BLambda = integral(ftheta, 0, pi);
        else
            % Calculate the normalization factor numerically
            Nsam = max(1000, 1000/P.DLambda);
            x1 = linspace(0, pi, Nsam);
            y1 = exp(-((sin(x1).^2 - P.Lambda0)./P.DLambda).^2).*sin(x1);
            BLambda1 = trapz(x1, y1);
            x2 = linspace(0, pi, 2*Nsam);
            y2 = exp(-((sin(x2).^2 - P.Lambda0)./P.DLambda).^2).*sin(x2);
            BLambda2 = trapz(x2, y2);
            if abs(BLambda1/BLambda2 - 1) <= errmax
                BLambda = BLambda2;
            else
                ftheta = @(x) exp(-((sin(x).^2 - P.Lambda0)./P.DLambda).^2).*sin(x);
                BLambda = integral(ftheta, 0, pi);
            end
        end
        Fac = 1.5./BLambda./log(1 + 1./uc.^3).*N.^4.*wp2/w2.*zeta.^2;
        sink = k_perp./ktot;
        sink2 = sink.^2;
        cosk = k_para./ktot;
        L0 = P.Lambda0;
        dL = P.DLambda;
        u = @(uyz) sqrt(zeta.^2 + uyz.^2);
        u2 = @(uyz) zeta.^2 + uyz.^2;
        up2 = @(uyz, theta) (zeta.*sink - uyz.*cosk.*cos(theta)).^2 + (uyz.*sin(theta)).^2;
        kup = @(uyz, theta) zeta.*sink2 - uyz.*sink.*cosk.*cos(theta);
        M = @(uyz, theta) exp(-((up2(uyz, theta)./u2(uyz) - L0)./dL).^2)./(u(uyz).^3 + uc.^3)...
            .*(3.*u(uyz).*zeta./(u(uyz).^3 + uc.^3) - 4.*(up2(uyz, theta)./u2(uyz) - L0)./(dL.^2)...
            .*(up2(uyz, theta).*zeta./(u2(uyz).^2) - kup(uyz, theta)./u2(uyz))).*uyz;
        ImD = integral2(M, 0, sqrt(1 - zeta.^2), 0, 2*pi, 'AbsTol', 1e-6);
        ImD = Fac.*ImD;
    end
end