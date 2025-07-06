function Pd = distp_xv(P, Ndist, Ps, EN, ifdf, Dim, kdir)
    % P is a struct given by loadparticle.m.
    % Ndist : number of sample points.
    % Ps : Particle struct array with fields vnorm(1*3).
    % Calculate distribution functions of each velocity element.

    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %         Note : This code is not carefully debugged!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if nargin < 7
        kdir = [1, 0, 0];
    end
    if nargin < 6
        Dim = [1, 2];
    end
    if nargin < 5
        ifdf = false;
    end
    if nargin < 4
        EN = 1;
    end
    if nargin < 3
        Ps.vnorm = [1, 1, 1];
        Ps = repmat(Ps, 1, max(1, numel(P.ele)) + max(1, numel(P.ion)));
    end
    if nargin < 2
        Ndist = 80;
    end
    Ndist = [Ndist(1), Ndist(end)];
    
    % Calculate eles
	Pe = Ps(1:EN);
    if ifdf
        numelPe = nan(size(P.ele));
        for k = 1:numel(P.ele)
            numelPe(k) = numel(P.ele(k).ft0)/sum(P.ele(k).ft0);
        end
    end
    for j = 1:numel(P.ele)
        if isempty(P.ele(j).v)
            continue;
        end
        Pd.ele(j).v = nan(3, Ndist(1));
        Pd.ele(j).r = nan(3, Ndist(end));
        Pd.ele(j).fxv = nan([Ndist, 3]);
        Pd.ele(j).fxv0 = nan([Ndist, 3]);
        Pd.ele(j).vperpnorm = nan;
        Pd.ele(j).vperp = nan(1, Ndist(1));
        Pd.ele(j).rperp = nan(1, Ndist(end));
        Pd.ele(j).fxv_perp = nan(Ndist);
        Pd.ele(j).fxv_perp0 = nan(Ndist);
        for k = 1:3
            if all(P.ele(j).r(k, :) == 0, 'all')
                X0 = kdir*P.ele(j).r;
            else
                X0 = P.ele(j).r(k, :);
            end
            if ifdf
                [Pd.ele(j).fxv(:, :, k), ve, xe] = histcountsdf2(P.ele(j).v(k, :), X0, P.ele(j).w, Ndist);
            	Pd.ele(j).v(k, :) = (ve(1:end-1) + ve(2:end))/2/Pe(j).vnorm(k);
				dx = (xe(end) - xe(1))/(numel(xe) - 1);
				Pd.ele(j).fxv(:, :, k) = Pd.ele(j).fxv(:, :, k)*numelPe(j)*Pe(j).vnorm(k);
                Pd.ele(j).fxv0(:, :, k) = repmat(exp(-(Pd.ele(j).v(k, :).').^2)*(1/sqrt(pi)/dx/Ndist(end)), 1, Ndist(end));
				Pd.ele(j).fxv(:, :, k) = Pd.ele(j).fxv(:, :, k) + Pd.ele(j).fxv0(:, :, k);
            else
                [Pd.ele(j).fxv(:, :, k), ve, xe] = histcounts2(P.ele(j).v(k, :), X0, Ndist, 'Normalization', 'pdf');
            	Pd.ele(j).v(k, :) = (ve(1:end-1) + ve(2:end))/2/Pe(j).vnorm(k);
				dx = (xe(end) - xe(1))/(numel(xe) - 1);
                Pd.ele(j).fxv(:, :, k) = Pd.ele(j).fxv(:, :, k)*Pe(j).vnorm(k);
                Pd.ele(j).fxv0(:, :, k) = repmat(exp(-(Pd.ele(j).v(k, :).').^2)*(1/sqrt(pi)/dx/Ndist(end)), 1, Ndist(end));
            end
            Pd.ele(j).r(k, :) = (xe(1:end-1) + xe(2:end))/2;
        end
        P.ele(j).vperp = sqrt(sum(P.ele.v(Dim, :).^2));
        Pd.ele(j).vperpnorm = 1./sqrt(sum(kdir./Pe(j).vnorm.^2));
        if ifdf
            [Pd.ele(j).fxv_perp, vperp, rperp] = histcountsdf2(P.ele(j).vperp, kdir*P.ele(j).r, P.ele(j).w, Ndist);
            Pd.ele(j).vperp = (vperp(1:end-1) + vperp(2:end))./2/Pd.ele(j).vperpnorm;
            Pd.ele(j).fxv_perp = Pd.ele(j).fxv_perp*Pd.ele(j).vperpnorm*numelPe(j);
            Pd.ele(j).fxv_perp0 = repmat(2/dx/Ndist(end)*(Pd.ele(j).vperp.*exp(-Pd.ele(j).vperp.^2)).', 1, Ndist(end));
            Pd.ele(j).fxv_perp = Pd.ele(j).fxv_perp + Pd.ele(j).fxv_perp0;
        else
            [Pd.ele(j).fxv_perp, vperp, rperp] = histcounts2(P.ele(j).vperp, kdir*P.ele(j).r, Ndist, 'Normalization', 'pdf');
            Pd.ele(j).vperp = (vperp(1:end-1) + vperp(2:end))./2/Pd.ele(j).vperpnorm;
            Pd.ele(j).fxv_perp = Pd.ele(j).fxv_perp*Pd.ele(j).vperpnorm;
            Pd.ele(j).fxv_perp0 = repmat(2/dx/Ndist(end)*(Pd.ele(j).vperp.*exp(-Pd.ele(j).vperp.^2)).', 1, Ndist(end));
        end
        Pd.ele(j).rperp = (rperp(1:end-1) + rperp(2:end))./2;
    end

    % Calculate ions
	Pi = Ps(EN+1:end);
    if ifdf
        numelPi = nan(size(P.ion));
        for k = 1:numel(P.ion)
            numelPi(k) = numel(P.ion(k).ft0)/sum(P.ion(k).ft0);
        end
    end
    for j = 1:numel(P.ion)
        if isempty(P.ion(j).v)
            continue;
        end
        Pd.ion(j).v = nan(3, Ndist(1));
        Pd.ion(j).r = nan(3, Ndist(end));
        Pd.ion(j).fxv = nan([Ndist, 3]);
        Pd.ion(j).fxv0 = nan([Ndist, 3]);
        Pd.ion(j).vperpnorm = nan;
        Pd.ion(j).vperp = nan(1, Ndist(1));
        Pd.ion(j).rperp = nan(1, Ndist(end));
        Pd.ion(j).fxv_perp = nan(Ndist);
        Pd.ion(j).fxv_perp0 = nan(Ndist);
        for k = 1:3
            if all(P.ion(j).r(k, :) == 0, 'all')
                X0 = kdir*P.ion(j).r;
            else
                X0 = P.ion(j).r(k, :);
            end
            if ifdf
                [Pd.ion(j).fxv(:, :, k), vi, xi] = histcountsdf2(P.ion(j).v(k, :), X0, P.ion(j).w, Ndist);
            	Pd.ion(j).v(k, :) = (vi(1:end-1) + vi(2:end))/2/Pi(j).vnorm(k);
				dx = (xe(end) - xe(1))/(numel(xe) - 1);
				Pd.ion(j).fxv(:, :, k) = Pd.ion(j).fxv(:, :, k)*numelPi(j)*Pi(j).vnorm(k);
                Pd.ion(j).fxv0(:, :, k) = repmat(exp(-(Pd.ion(j).v(k, :).').^2)*(1/sqrt(pi)/dx/Ndist(end)), 1, Ndist(end));
				Pd.ion(j).fxv(:, :, k) = Pd.ion(j).fxv(:, :, k) + Pd.ion(j).fxv0(:, :, k);
            else
                [Pd.ion(j).fxv(:, :, k), vi, xi] = histcounts2(P.ion(j).v(k, :), X0, Ndist, 'Normalization', 'pdf');
                Pd.ion(j).v(k, :) = (vi(1:end-1) + vi(2:end))/2/Pi(j).vnorm(k);
				dx = (xe(end) - xe(1))/(numel(xe) - 1);
                Pd.ion(j).fxv(:, :, k) = Pd.ion(j).fxv(:, :, k)*Pi(j).vnorm(k);
                Pd.ion(j).fxv0(:, :, k) = repmat(exp(-(Pd.ion(j).v(k, :).').^2)*(1/sqrt(pi)/dx/Ndist(end)), 1, Ndist(end));
            end
            Pd.ion(j).r(k, :) = (xi(1:end-1) + xi(2:end))/2;
        end
        P.ion(j).vperp = sqrt(sum(P.ion(j).v(Dim, :).^2));
        Pd.ion(j).vperpnorm = 1./sqrt(sum(kdir./Pi(j).vnorm.^2));
        if ifdf
            [Pd.ion(j).fxv_perp, vperp, rperp] = histcountsdf2(P.ion(j).vperp, kdir*P.ion(j).r, P.ion(j).w, Ndist);
            Pd.ion(j).vperp = (vperp(1:end-1) + vperp(2:end))./2/Pd.ion(j).vperpnorm;
            Pd.ion(j).fxv_perp = Pd.ion(j).fxv_perp*Pd.ion(j).vperpnorm*numelPi(j);
            Pd.ion(j).fxv_perp0 = repmat(2/dx/Ndist(end)*(Pd.ion(j).vperp.*exp(-Pd.ion(j).vperp.^2)).', 1, Ndist(end));
            Pd.ion(j).fxv_perp = Pd.ion(j).fxv_perp + Pd.ion(j).fxv_perp0;
        else
            [Pd.ion(j).fxv_perp, vperp, rperp] = histcounts2(P.ion(j).vperp, kdir*P.ion(j).r, Ndist, 'Normalization', 'pdf');
            Pd.ion(j).vperp = (vperp(1:end-1) + vperp(2:end))./2/Pd.ion(j).vperpnorm;
            Pd.ion(j).fxv_perp = Pd.ion(j).fxv_perp*Pd.ion(j).vperpnorm;
            Pd.ion(j).fxv_perp0 = repmat(2/dx/Ndist(end)*(Pd.ion(j).vperp.*exp(-Pd.ion(j).vperp.^2)).', 1, Ndist(end));
        end
        Pd.ion(j).rperp = (rperp(1:end-1) + rperp(2:end))./2;
    end
end
