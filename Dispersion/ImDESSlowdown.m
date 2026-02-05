function ImD = ImDESSlowdown(w, k_para, k_perp, B, P, ifautomatic)
% Calculate Imaginary Part of Slowing-down distribution.
% ImD = N^2 N.\Chi^a.N
    errmax = 5e-3;
    if nargin < 6
        ifautomatic = false;
    end
    % -------------------- Physical Constants -------------------- %
    c = 299792458;
    e = 1.602176565e-19;
    epsilon_0 = 8.854187817e-12;

    % -------------------- Normalization Factors -------------------- %
    % Calculate the normalization factor with fixed step
    isotropic = ~isfield(P, 'DLambda') || P.DLambda == 0 || isnan(P.DLambda) || isempty(P.DLambda);
    if isotropic
        BLambda = 2;
    else
        if ifautomatic
            % Calculate the normalization factor automatically
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
    end

    % -------------------- Some Factors -------------------- %
    ktot2 = k_para^2 + k_perp^2;
    N2 = ktot2*(c/w)^2;
    wc = P.q.*e.*B./P.m;
    wp2 = (P.q.*e).^2.*P.n0./epsilon_0./P.m;
    uc = P.uc./P.u0;
    uc3 = uc.^3;
    rhoperp = k_perp.*P.u0./wc;
    rhopara = k_para.*P.u0./wc;
    zeta_perp = w./k_perp./P.u0;
    zeta_para = abs(w./k_para./P.u0);
    % Factor for Im(D)
    Fac = (wp2./w.^2).*zeta_perp.^2.*zeta_para.*(9*pi./log(1 + (P.u0./P.uc).^3)./BLambda);

    % -------------------- Calculate Integral -------------------- %
    nmin = ceil((w - abs(k_para).*P.u0)./abs(wc));
    nmax = floor((w + abs(k_para).*P.u0)./abs(wc));
    IndN = nmin:nmax;
    un_para = (1./k_para./P.u0).*(w - wc*IndN);
    up_max = sqrt(1 - un_para.^2);
    if isotropic
        D = nan(size(un_para));
        if ifautomatic
            for n = 1:numel(D)
                % Calculate the intrgral automatically
                df = @(x) x.*sqrt(x.^2+un_para(n).^2)./((uc3 + (x.^2+un_para(n).^2).^1.5).^2).*besselj(IndN(n), rhoperp.*x).^2;
                D(n) = integral(df, 0, up_max(n), 'AbsTol', 1e-6);
            end
        else
            for n = 1:numel(D)
                % Calculate the intrgral with fixed step
                dx = 0.025*min([uc, 1./rhoperp]);
                Nsam = max(ceil(up_max(n)./dx), 1000);
                x1 = linspace(0, up_max(n), Nsam);
                x2 = linspace(0, up_max(n), 2*Nsam);
                % Calculate D
                y1 = x1.*sqrt(x1.^2+un_para(n).^2)./((uc3 + (x1.^2+un_para(n).^2).^1.5).^2).*besselj(IndN(n), rhoperp.*x1).^2;
                D1 = trapz(x1, y1);
                y2 = x2.*sqrt(x2.^2+un_para(n).^2)./((uc3 + (x2.^2+un_para(n).^2).^1.5).^2).*besselj(IndN(n), rhoperp.*x2).^2;
                D2 = trapz(x2, y2);
                if abs(D1/D2 - 1) < errmax
                    D(n) = D2;
                else
                    df = @(x) x.*sqrt(x.^2+un_para(n).^2)./((uc3 + (x.^2+un_para(n).^2).^1.5).^2).*besselj(IndN(n), rhoperp.*x).^2;
                    D(n) = integral(df, 0, up_max(n), 'AbsTol', 1e-6);
                end
            end
        end
        ImD = Fac.*sum(D, "all");
    else
        D = nan(size(un_para));
        E = nan(size(un_para));
        if ifautomatic
            for n = 1:numel(D)
                % Calculate the intrgral automatically
                df = @(x) x.*sqrt(x.^2+un_para(n).^2)./((uc3 + (x.^2+un_para(n).^2).^1.5).^2).*...
                    exp(-((x.^2./(x.^2+un_para(n).^2) - P.Lambda0)./P.DLambda).^2).*besselj(IndN(n), rhoperp.*x).^2;
                D(n) = integral(df, 0, up_max(n), 'AbsTol', 1e-6);
                ef = @(x) x.*(x.^2./(x.^2+un_para(n).^2) - P.Lambda0)./(P.DLambda.^2)./(uc3 + (x.^2+un_para(n).^2).^1.5)./((x.^2+un_para(n).^2).^2).*...
                    un_para(n).*(-IndN(n).*un_para(n) + rhopara.*x.^2).*exp(-((x.^2./(x.^2+un_para(n).^2) - P.Lambda0)./P.DLambda).^2).*besselj(IndN(n), rhoperp.*x).^2;
                E(n) = integral(ef, 0, up_max(n), 'AbsTol', 1e-6);
            end
        else
            for n = 1:numel(D)
                % Calculate the intrgral with fixed step
                dx = 0.025*min([uc, P.DLambda./2, 1./rhoperp]);
                Nsam = max(ceil(up_max(n)./dx), 1000);
                x1 = linspace(0, up_max(n), Nsam);
                x2 = linspace(0, up_max(n), 2*Nsam);
                % Calculate D
                y1 = x1.*sqrt(x1.^2+un_para(n).^2)./((uc3 + (x1.^2+un_para(n).^2).^1.5).^2).*...
                    exp(-((x1.^2./(x1.^2+un_para(n).^2) - P.Lambda0)./P.DLambda).^2).*besselj(IndN(n), rhoperp.*x1).^2;
                D1 = trapz(x1, y1);
                y2 = x2.*sqrt(x2.^2+un_para(n).^2)./((uc3 + (x2.^2+un_para(n).^2).^1.5).^2).*...
                    exp(-((x2.^2./(x2.^2+un_para(n).^2) - P.Lambda0)./P.DLambda).^2).*besselj(IndN(n), rhoperp.*x2).^2;
                D2 = trapz(x2, y2);
                if abs(D1/D2 - 1) < errmax
                    D(n) = D2;
                else
                    df = @(x) x.*sqrt(x.^2+un_para(n).^2)./((P.uc.^3 + (x.^2+un_para(n).^2).^1.5).^2).*...
                        exp(-((x.^2./(x.^2+un_para(n).^2) - P.Lambda0)./P.DLambda).^2).*besselj(IndN(n), rhoperp.*x).^2;
                    D(n) = integral(df, 0, up_max(n), 'AbsTol', 1e-6);
                end
                % Calculate E
                z1 = x1.*(x1.^2./(x1.^2+un_para(n).^2) - P.Lambda0)./(P.DLambda.^2)./(uc3 + (x1.^2+un_para(n).^2).^1.5)./((x1.^2 + un_para(n).^2).^2).*...
                    un_para(n).*(-IndN(n).*un_para(n) + rhopara.*x1.^2).*exp(-((x1.^2./(x1.^2+un_para(n).^2) - P.Lambda0)./P.DLambda).^2).*besselj(IndN(n), rhoperp.*x1).^2;
                E1 = trapz(x1, z1);
                z2 = x2.*(x2.^2./(x2.^2+un_para(n).^2) - P.Lambda0)./(P.DLambda.^2)./(uc3 + (x2.^2+un_para(n).^2).^1.5)./((x2.^2 + un_para(n).^2).^2).*...
                    un_para(n).*(-IndN(n).*un_para(n) + rhopara.*x2.^2).*exp(-((x2.^2./(x2.^2+un_para(n).^2) - P.Lambda0)./P.DLambda).^2).*besselj(IndN(n), rhoperp.*x2).^2;
                E2 = trapz(x2, z2);
                if abs(E1./E2 - 1) < errmax
                    E(n) = E2;
                else
                    ef = @(x) x.*(x.^2./(x.^2+un_para(n).^2) - P.Lambda0)./(P.DLambda.^2)./(uc3 + (x.^2+un_para(n).^2).^1.5)./((x.^2 + un_para(n).^2).^2).*...
                        un_para(n).*(-IndN(n).*un_para(n) + rhopara.*x.^2).*exp(-((x.^2./(x.^2+un_para(n).^2) - P.Lambda0)./P.DLambda).^2).*besselj(IndN(n), rhoperp.*x).^2;
                    E(n) = integral(ef, 0, up_max(n), 'AbsTol', 1e-6);
                end
            end
        end
        ImD = Fac.*(sum(D, "all") - (4*wc./3./w).*sum(E, "all"));
    end
    ImD = N2.*(k_perp*c./w).^2.*ImD;
end