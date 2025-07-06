function varargout = DBram1976(w, k_para, k_perp, B, Ps, EN)
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
    N_para = k_para*c./w;
    N_perp = k_perp*c./w;
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
    Dis = zeros(1, numel(Ps) - EN);
    Dis0 = zeros(numel(Ps) - EN, 1);
    Disp = zeros(numel(Ps) - EN, 500);
    Disn = zeros(numel(Ps) - EN, 500);
    Disa = zeros(numel(Ps) - EN, 1001);
    for k = EN+1 : numel(Ps)
        Dik1 = 1 + Ps(k).Lperp.*real(Z_plasma(Ps(k).Lperp));
        % Calculate Summary
        %Summation start at the gyro hamonic
        ngyro = round(w./Ps(k).wc);
        Y = Ps(k).wc./k_para./Ps(k).vt;
        Dik2 = dZ_plasma(Ps(k).Lpara - ngyro.*Y)./(w - ngyro.*Ps(k).wc);
        Dis0(k-EN) = Dik2;
        n = 1;
        Disp(k-EN, n) = dZ_plasma(Ps(k).Lpara - (ngyro + n).*Y)./(w - (ngyro + n).*Ps(k).wc);
        Disn(k-EN, n) = dZ_plasma(Ps(k).Lpara - (ngyro - n).*Y)./(w - (ngyro - n).*Ps(k).wc);
        dDik2 = Disp(k-EN, n) + Disn(k-EN, n);
        % dDik2 = dZ_plasma(Ps(k).Lpara - (ngyro + n).*Y)./(w - (ngyro + n).*Ps(k).wc)...
        %         + dZ_plasma(Ps(k).Lpara - (ngyro - n).*Y)./(w - (ngyro - n).*Ps(k).wc);
        ifbreak = abs(real(dDik2)) <= epse*abs(real(Dik2)) & abs(imag(dDik2)) <= epse*abs(imag(Dik2));
        % ifbreak = abs(dDik2) <= epse*abs(Dik2);
        while(~ifbreak)
            n = n + 1;
            Dik2 = Dik2 + dDik2;
            Disp(k-EN, n) = dZ_plasma(Ps(k).Lpara - (ngyro + n).*Y)./(w - (ngyro + n).*Ps(k).wc);
            Disn(k-EN, n) = dZ_plasma(Ps(k).Lpara - (ngyro - n).*Y)./(w - (ngyro - n).*Ps(k).wc);
            dDik2 = Disp(k-EN, n) + Disn(k-EN, n);
            % dDik2 = dZ_plasma(Ps(k).Lpara - (ngyro + n).*Y)./(w - (ngyro + n).*Ps(k).wc)...
            %         + dZ_plasma(Ps(k).Lpara - (ngyro - n).*Y)./(w - (ngyro - n).*Ps(k).wc);
            ifbreak = abs(real(dDik2)) <= epse*abs(real(Dik2)) & abs(imag(dDik2)) <= epse*abs(imag(Dik2));
            % ifbreak = abs(dDik2) <= epse*abs(Dik2);
            if n > 500 %10000
                fprintf('Summary fail in Brambilla\n');
                D = NaN;
                break;
            end
        end
        Disa(k-EN, :) = [fliplr(Disn(k-EN, :)), Dis0(k-EN), Disp(k-EN, :)];
        Disa(k-EN, :) = Ps(k).wc./2./pi.*Disa(k-EN, :);
        Disa(k-EN, :) = -1.772453850905516.*Ps(k).Lperp.*exp(-Ps(k).Lperp.^2).*Disa(k-EN, :);
        Disa(k-EN, :) = 2*(Ps(k).wp.*c./w./Ps(k).vt).^2.*Disa(k-EN, :);

        Dik2 = Dik2 + dDik2;
        Dik2 = cot(pi*w./Ps(k).wc) + Ps(k).wc./2./pi.*Dik2;
        Dik2 = -1.772453850905516.*Ps(k).Lperp.*exp(-Ps(k).Lperp.^2).*Dik2;
        Dik = 2*(Ps(k).wp.*c./w./Ps(k).vt).^2.*(Dik1 + Dik2);
        Di = Di + Dik;
        Dis(k - EN) = Dik;
    end
    D = Di + De;
    Ds = [Des, Dis];
    varargout{1} = D;
    varargout{2} = Ds.';
    varargout{3} = Disa;
    %--------------Check Units--------------%
    N2 = N_para.^2 + N_perp.^2;
    varargout{1} = varargout{1}.*N2;
    varargout{2} = varargout{2}.*N2;
    varargout{3} = Disa.*N2;
end

%% Special functions
%--------------Plasma dispersion function--------------%
function dZ = dZ_plasma(z)
    dZ = -2.*(1 + z.*Z_plasma(z));
end