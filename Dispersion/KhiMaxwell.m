function Khi = KhiMaxwell(w, k_para, k_perp, B, P, Psi)
% Calculate susceptibility tensor for single P species.
% MATLAB functions Z_plasma.m and cef.m are required.
    if nargin < 6
        Psi = 0;
    end
    %--------------Default Input parameters--------------%
    errmax = 1e-7;
    errmaxI = 5e-7;        %Relative error compared to unit.
    maxsumnum = 10000;
    %--------------Physical constants--------------%
    epsilon_0 = 8.854187817e-12;
    e = 1.602176565e-19;
    %--------------Add Temperature Unisotropicity and Drift--------------%
    if ~isfield(P, 'T_para')
        P.T_para = P.T;
    else
        if isempty(P.T_para)
            P.T_para = P.T;
        end
    end
    if ~isfield(P, 'T_perp')
        P.T_perp = P.T;
    else
        if isempty(P.T_perp)
            P.T_perp = P.T;
        end
    end
    if ~isfield(P, 'v0')
        P.v0 = 0;
    else
        if isempty(P.v0)
            P.v0 = 0;
        end
    end
    %--------------Frequently used parameters--------------%
    % == Group 1
    sig = sign(P.q);
    wc = sig.*P.q.*e.*B./P.m;%Always postive.
    wp2 = (P.q.*e).^2.*P.n0./(epsilon_0.*P.m);
    v_para = sqrt(2.*P.T_para.*e./P.m);
    v_perp = sqrt(2.*P.T_perp.*e./P.m);
    rho_L = v_perp./wc;
    % == Group 2
    X = wp2./(w.^2);
    r_para = k_para.*v_para./w;
    r_0 = k_para.*P.v0./w;
    lambda = 0.5.*(k_perp.*rho_L).^2;
    Y = wc./w;
    mu = P.T_perp./P.T_para;
    nu = k_perp./k_para;
    %--------------Cofficients used after the summary--------------%
    Cff = X.*exp(-1i*imag(lambda)).*[(1./r_para).*[2.*lambda;     1./lambda;     1i.*sig;     -1];       (nu./Y).*[1./lambda;     1i.*sig]];
    %--------------Summary over n--------------%
    % == Calculate gyro resonance
    %Summation will be performed from 0 and ngyro simutaneously
    ngyro = - round(1./abs(Y));  %Note that wc > 0; zeta = (1 + n.*Y)./r_para
    maxsumnum1 = max(maxsumnum, 20*abs(ngyro));
    % == Initialize
    N = 0;                                              %Number of terms added from both points.
    dN = max(10, round(0.1*abs(ngyro)));                %Add dN terms every step.
    k_zero = CalcK(N, r_para, r_0, lambda, Y, mu);
    if ngyro ~= 0
        k_gyro = CalcK(ngyro + N, r_para, r_0, lambda, Y, mu);
    else
        k_gyro = zeros(size(k_zero, 1), 1);
    end
    k = sum([k_zero, k_gyro], 2);
    % == Calculate increments
    if ngyro + N + dN < 0 - N - dN
        Nind = [-N-dN:-N-1, N+1:N+dN];
        k_zero = CalcK(0 + Nind, r_para, r_0, lambda, Y, mu);
        k_gyro = CalcK(ngyro + Nind, r_para, r_0, lambda, Y, mu);
        dk = sum([k_zero, k_gyro], 2);
    else
        Nindl = -N-dN:-N-1;
        Nindr = N+1:N+dN;
        Nindm = ngyro+N+1 : 0-N-1;
        if ~isempty(Nindm)
            k_mid = CalcK(Nindm, r_para, r_0, lambda, Y, mu);
        else
            k_mid = zeros(6, 1);
        end
        k_zero = CalcK(0 + Nindr, r_para, r_0, lambda, Y, mu);
        k_gyro = CalcK(ngyro + Nindl, r_para, r_0, lambda, Y, mu);
        dk = sum([k_zero, k_mid, k_gyro], 2);
    end
    % == Check error
    errmax1 = errmaxI/Cff(2);
    errmax3 = errmaxI/Cff(4);
    %Array to show whether elements in k meet the error limit.
    iserrR = abs(real(dk)) <= errmax*abs(real(k));
    iserrI = abs(imag(dk)) <= errmax*abs(imag(k));
    iserr = iserrR & iserrI;
    %K1 and K3 are compared to 1
    if abs(k(2)) <= errmax1
        iserr(2) = 1;
    end
    if abs(k(4)) <= errmax3
        iserr(4) = 1;
    end
    %Cases : 0 <= 0*errmax   -----   true, meets the error limit
    %        1 <= 0*errmax   -----   false, exceeds the error limit
    %        0 <= 1*errmax   -----   true, meets the error limit
    %        1 <= 1*errmax   -----   false, exceeds the error limit
    %Signal showing whether the summation shoule continue.
    %Continue if not all of the elements meet the error limit.
    ifcontinue = sum(iserr, 'all') < numel(iserr);
    % == Perform summation in a loop
    while(ifcontinue)
        k = k + dk;
        N = N + dN;
        % == Calculate increments
        if ngyro + N + dN < 0 - N - dN
            Nind = [-N-dN:-N-1, N+1:N+dN];
            k_zero = CalcK(0 + Nind, r_para, r_0, lambda, Y, mu);
            k_gyro = CalcK(ngyro + Nind, r_para, r_0, lambda, Y, mu);
            dk = sum([k_zero, k_gyro], 2);
        else
            Nindl = -N-dN:-N-1;
            Nindr = N+1:N+dN;
            Nindm = ngyro+N+1 : 0-N-1;
            if ~isempty(Nindm)
                k_mid = CalcK(Nindm, r_para, r_0, lambda, Y, mu);
            else
                k_mid = zeros(6, 1);
            end
            k_zero = CalcK(0 + Nindr, r_para, r_0, lambda, Y, mu);
            k_gyro = CalcK(ngyro + Nindl, r_para, r_0, lambda, Y, mu);
            dk = sum([k_zero, k_mid, k_gyro], 2);
        end
        % == Check error
        %Array to show whether elements in k meet the error limit.
        iserrR = abs(real(dk)) <= errmax*abs(real(k));
        iserrI = abs(imag(dk)) <= errmax*abs(imag(k));
        iserr = iserrR & iserrI;
        %K1 and K3 are compared to 1
        if abs(k(2)) <= errmax1
            iserr(2) = 1;
        end
        if abs(k(4)) <= errmax3
            iserr(4) = 1;
        end
        %Cases : 0 <= 0*errmax   -----   true, meets the error limit
        %        1 <= 0*errmax   -----   false, exceeds the error limit
        %        0 <= 1*errmax   -----   true, meets the error limit
        %        1 <= 1*errmax   -----   false, exceeds the error limit
        %Signal showing whether the summation shoule continue.
        %Continue if not all of the elements meet the error limit.
        ifcontinue = sum(iserr, 'all') < numel(iserr);
        if N > maxsumnum1 && ifcontinue
            % break;
            % fprintf('--------------------------Summation diverges in KhiMaxwell.m--------------------------\n');
            % fprintf('Terms added : %d\n', N);
            % fprintf('RelErr : %f\n', errmax);
            % fprintf('Convergence : \n');
            % disp(abs(dk).');
            % disp(abs(k).');
            % disp(iserr.');
            % fprintf('w, k_para, k_perp, B:\n');
            % disp([w, k_para, k_perp, B]);
            % fprintf('--------------------------End of error massage in KhiMaxwell.m--------------------------\n');
            Khi = nan(3);
            return;
        end
    end
    K = k + dk;
    %--------------Calculate K0~K5--------------%
    K = Cff.*K;
    %--------------Derive the matrix and rotate--------------%
	Khi = [K(2),    K(3),           K(5);...
             -K(3),   K(2)+K(1),  -K(6);...
             K(5),    K(6),           K(4)];
    Cos = cos(Psi);
    Sin = sin(Psi);
    T = [Cos,   -Sin,   0;...
           Sin,    Cos,   0;...
           0,       0,      1];
    Khi = T*Khi*(T.');
end

function K = CalcK(n, r_para, r_0, lambda, Y, mu)
% Calculate single terms in the summary of k0~k5.
% n is a row vector.
    zeta = (1 - r_0 + n.*Y)./r_para;
    %Plasma dispersion function and its derivative
    Z_zeta = Z_plasma(zeta);
    Z_zeta_1 = -2.*(1 + zeta.*Z_zeta);
    %First modified Bessel function and its derivative
    In = besseli(n, lambda, 1);
    In1 = (besseli(n-1, lambda, 1) + besseli(n+1, lambda, 1))./2;
    %Common usage of plasma dispersion function
    k_02 = (1 - r_0).*Z_zeta + 0.5.*r_para.*(1 - mu).*Z_zeta_1;
    if r_0 ~= 0
        k_3 = ((1 + n.*Y)./r_para).*((1 + n.*Y.*(1 - 1./mu)).*Z_zeta_1 + 2.*n.*(Y./mu).*(r_0./r_para).*(Z_zeta + r_para./(1 + n.*Y)));
    else
        k_3 = ((1 + n.*Y)./r_para).*((1 + n.*Y.*(1 - 1./mu)).*Z_zeta_1);
    end
    k_45 = n.*Y.*(r_0./r_para).*Z_zeta + 0.5.*(mu - n.*Y.*(1 - mu)).*Z_zeta_1;
    %k0~k5
    k0 = (In - In1).*k_02;
    k1 = n.^2.*In.*k_02;
    k2 = n.*(In - In1).*k_02;
    k3 = In.*k_3;
    k4 = n.*In.*k_45;
    k5 = (In - In1).*k_45;
    K = [k0; k1; k2; k3; k4; k5];
end
