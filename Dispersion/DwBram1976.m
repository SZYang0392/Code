function varargout = DwBram1976(w, k_para, k_perp, B, Ps, EN)
% Calculate Brambilla Dispersion
% w, N_para, N_perp, B, Ps(k).* can be vectors of the same size.
% Ps(1:EN) denote electrons. The remaining denote ions.
    if nargin < 6
        EN = 1;
    end
    %--------------Constants list--------------%
    e = 1.602176565e-19;
    c = 299792458;
    % e^2/epsi0 = 2.899158905e-27
    % sqrt(2e) = 5.660700602e-10
    % sqrt(pi) = 1.772453850905516
    %--------------Refractive Index--------------%
    N_para = k_para.*c./w;
    N_perp = k_perp.*c./w;
    %--------------Particle Parameters--------------%
    for k = 1:numel(Ps)
        Ps(k).wp2 = 2.899158905e-27.*Ps(k).q.^2.*Ps(k).n0./Ps(k).m;
        Ps(k).wp = sqrt(Ps(k).wp2);
        Ps(k).wc = e.*abs(Ps(k).q).*B./Ps(k).m;
        Ps(k).wc2 = Ps(k).wc.^2;
        Ps(k).vt = 5.660700602e-10*sqrt(Ps(k).T./Ps(k).m);
        Ps(k).Lpara = w./k_para./Ps(k).vt;
        Ps(k).Lperp = w./k_perp./Ps(k).vt;
    end

    %%===================Calculate D===================%%
    %--------------Electrons--------------%
    Des = zeros(1, EN);
    for k = 1:EN
        Des(k) = Des(k) + N_perp.^2.*(1 + Ps(k).wp2./Ps(k).wc2);
        Des(k) = Des(k) - N_para.^2.*(Ps(k).wp./w.*Ps(k).Lpara).^2.*dZ_plasma(Ps(k).Lpara);
    end
    De = sum(Des, 'all');
    %--------------Ions--------------%
    epse = 1e-5;
    Di = 0;
    Dik1s = zeros(1, numel(Ps) - EN);
    Dik2s = zeros(1, numel(Ps) - EN);
    ns = zeros(1, numel(Ps) - EN);
    Dis = zeros(1, numel(Ps) - EN);
    for k = EN+1 : numel(Ps)
        Dik1 = 1 + Ps(k).Lperp.*real(Z_plasma(Ps(k).Lperp));
        Dik1s(k-EN) = Dik1;
        % Calculate Summary
        %Summation start at the gyro hamonic
        ngyro = round(w./Ps(k).wc);
        Y = Ps(k).wc./k_para./Ps(k).vt;
        Dik2 = dZ_plasma(Ps(k).Lpara - ngyro.*Y)./(w - ngyro.*Ps(k).wc);
        n = 1;
        dDik2 = dZ_plasma(Ps(k).Lpara - (ngyro + n).*Y)./(w - (ngyro + n).*Ps(k).wc)...
                + dZ_plasma(Ps(k).Lpara - (ngyro - n).*Y)./(w - (ngyro - n).*Ps(k).wc);
        ifbreak = abs(real(dDik2)) <= epse*abs(real(Dik2)) & abs(imag(dDik2)) <= epse*abs(imag(Dik2));
        % ifbreak = abs(dDik2) <= epse*abs(Dik2);
        while(~ifbreak)
            n = n + 1;
            Dik2 = Dik2 + dDik2;
            dDik2 = dZ_plasma(Ps(k).Lpara - (ngyro + n).*Y)./(w - (ngyro + n).*Ps(k).wc)...
                    + dZ_plasma(Ps(k).Lpara - (ngyro - n).*Y)./(w - (ngyro - n).*Ps(k).wc);
            ifbreak = abs(real(dDik2)) <= epse*abs(real(Dik2)) & abs(imag(dDik2)) <= epse*abs(imag(Dik2));
            % ifbreak = abs(dDik2) <= epse*abs(Dik2);
            if n > 7500 & ~ifbreak
                fprintf('Summary fail in Brambilla\n');
                D = NaN;
                break;
            end
        end
        ns(k) = n;
        Dik2 = Dik2 + dDik2;
        Dik2 = cot(pi*w./Ps(k).wc) + Ps(k).wc./2./pi.*Dik2;
        Dik2 = -1.772453850905516.*Ps(k).Lperp.*exp(-Ps(k).Lperp.^2).*Dik2;
        Dik2s(k-EN) = Dik2;
        Dik = 2*(Ps(k).wp.*c./w./Ps(k).vt).^2.*(Dik1 + Dik2);
        Di = Di + Dik;
        Dis(k-EN) = Dik;
    end
    D = Di + De;
    varargout{1} = D;
    varargout{3} = [Des, Dis].';
    %%===================Calculate DwD===================%%
    % Note : the summation has an identical truncation to that in calculating D.
    DwD1 = -2*D./w;
    %--------------Electrons--------------%
    DwD2 = 0;
    for k = 1 : EN
        DwD2 = DwD2 - Ps(k).wp2./(Ps(k).vt.^3).*dZ_plasma2(Ps(k).Lpara);
    end
    DwD2 = DwD2.*N_para.^2./k_para.^3;
    DwDe = DwD2;
    %--------------Ions--------------%
    DwD3 = 0;
    DwD4 = 0;
    DwD5 = 0;
    DwD6 = 0;
    DwD7 = 0;
    DwD8 = 0;
    for k = EN+1 : numel(Ps)
        param = 2.*(Ps(k).wp.*c./w./Ps(k).vt).^2;
        DwD3 = DwD3 + param.*(Dik1s(k-EN) + Dik2s(k-EN) - 1)./w;
        DwD4 = DwD4 + param.*(-2./w).*Ps(k).Lperp.^2.*Dik2s(k-EN);
        param1 = param.*1.772453850905516.*Ps(k).Lperp.*exp(-Ps(k).Lperp.^2);
        DwD5 = DwD5 - param1 .*(-1)./(sin(pi*w./Ps(k).wc).^2).*pi./Ps(k).wc;
        ngyro = round(w./Ps(k).wc);
        n = -ns(k-EN):ns(k-EN) + ngyro;
        Y = Ps(k).wc./k_para./Ps(k).vt;
        zeta = Ps(k).Lpara - n.*Y;
        DwD60 = -sum(dZ_plasma(zeta)./((w - n.*Ps(k).wc).^2), 2);
        DwD6 = DwD6 - param1  .*Ps(k).wc./2./pi.*DwD60;
        DwD70 = sum(dZ_plasma2(zeta)./(w - n.*Ps(k).wc), 2);
        DwD7 = DwD7 - param1  .*Ps(k).Lpara./w.*Ps(k).wc./2./pi.*DwD70;
        DwD8 = DwD8 + param.*Ps(k).Lperp.^2./w.*real(dZ_plasma(Ps(k).Lperp));
    end
    DwDi = DwD3 + DwD4 + DwD5 + DwD6 + DwD7 + DwD8;
    DwD = DwD1 + DwDe + DwDi;
    varargout{2} = DwD;
    %--------------Check Units--------------%
    N2 = N_para.^2 + N_perp.^2;
    varargout{1} = varargout{1}.*N2;
    varargout{2} = varargout{2}.*N2 + (-2./w).*varargout{1};
    varargout{3} = varargout{3}.*N2;
end

%% Special functions
%--------------Plasma dispersion function--------------%
function dZ = dZ_plasma2(z)
    dZ = -2.*(Z_plasma(z) + z.*dZ_plasma(z));
end
function dZ = dZ_plasma(z)
    dZ = -2.*(1 + z.*Z_plasma(z));
end