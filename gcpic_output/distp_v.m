function Pd = distp_v(P, Ndist, Ps, EN, ifdf, Dir)
    % P is a struct given by loadparticle.m.
    % Ndist : number of sample points.
    % Ps : Particle struct array with fields vnorm(1*3).
    % Calculate distribution functions of each velocity element.

    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %         Note : This code is not carefully debugged!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if nargin >= 6
        for j = 1:numel(P.ele)
            P.ele(j).r = (Dir.')*P.ele(j).r;
            if isfield(P.ele(j), vb)
                P.ele(j).v = P.ele(j).vb;
            else
                P.ele(j).v = (Dir.')*P.ele(j).v;
            end
        end
        for j = 1:numel(P.ion)
            P.ion(j).r = (Dir.')*P.ion(j).r;
            if isfield(P.ion(j), vb)
                P.ion(j).v = P.ion(j).vb;
            else
                P.ion(j).v = (Dir.')*P.ion(j).v;
            end
        end
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
    vsam = nan(3, Ndist + 1);
    
    % Calculate eles
    if nargin >= 3
		Pe = Ps(1:EN);
    end
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
        Pd.ele(j).v = nan(3, Ndist);
        Pd.ele(j).fv = nan(3, Ndist);
        Pd.ele(j).fv0 = nan(3, Ndist);
        for k = 1:3
            if ifdf
                [Pd.ele(j).fv(k, :), vsam(k, :)] = histcountsdf(P.ele(j).v(k, :), P.ele(j).w, Ndist);
				Pd.ele(j).v(k, :) = (vsam(k, 1:end-1) + vsam(k, 2:end))/2/Pe(j).vnorm(k);
                Pd.ele(j).fv(k, :) = Pd.ele(j).fv(k, :)*numelPe(j)*Pe(j).vnorm(k);
                Pd.ele(j).fv0(k, :) = exp(-Pd.ele(j).v(k, :).^2)/sqrt(pi);
				Pd.ele(j).fv(k, :) = Pd.ele(j).fv(k, :) + Pd.ele(j).fv0(k, :);
            else
                [Pd.ele(j).fv(k, :), vsam(k, :)] = histcounts(P.ele(j).v(k, :), Ndist, 'Normalization', 'pdf');
				Pd.ele(j).v(k, :) = (vsam(k, 1:end-1) + vsam(k, 2:end))/2;
                Pd.ele(j).fv0(k, :) = exp(-Pd.ele(j).v(k, :).^2)/sqrt(pi);
		        if nargin >= 3
		            Pd.ele(j).fv(k, :) = Pd.ele(j).fv(k, :)*Pe(j).vnorm(k);
		            Pd.ele(j).v(k, :) = Pd.ele(j).v(k, :)/Pe(j).vnorm(k);
		        end
            end
        end
        if ifdf
            [Pd.ele(j).fvperp, ~, ~] = histcountsdf2(P.ele(j).v(2, :), P.ele(j).v(1, :), P.ele(j).w, [Ndist, Ndist], vsam(2, :), vsam(1, :));
            Pd.ele(j).fvperp = Pd.ele(j).fvperp*Pe(j).vnorm(1)*Pe(j).vnorm(2)*numelPe(j);
            Pd.ele(j).fvperp0 = exp(-(Pd.ele(j).v(2, :).').^2 - Pd.ele(j).v(1, :).^2)/pi;
            Pd.ele(j).fvperp = Pd.ele(j).fvperp + Pd.ele(j).fvperp0;
        else
            [Pd.ele(j).fvperp, ~, ~] = histcounts2(P.ele(j).v(2, :), P.ele(j).v(1, :), vsam(2, :), vsam(1, :), 'Normalization', 'pdf');
            Pd.ele(j).fvperp = Pd.ele(j).fvperp*Pe(j).vnorm(1)*Pe(j).vnorm(2);
            Pd.ele(j).fvperp0 = exp(-(Pd.ele(j).v(2, :).').^2 - Pd.ele(j).v(1, :).^2)/pi;
        end
    end

    % Calculate ions
    if nargin >= 3
		Pi = Ps(EN+1:end);
    end
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
        Pd.ion(j).v = nan(3, Ndist);
        Pd.ion(j).fv = nan(3, Ndist);
        Pd.ion(j).fv0 = nan(3, Ndist);
        for k = 1:3
            if ifdf
                [Pd.ion(j).fv(k, :), vsam(k, :)] = histcountsdf(P.ion(j).v(k, :), P.ion(j).w, Ndist);
				Pd.ion(j).v(k, :) = (vsam(k, 1:end-1) + vsam(k, 2:end))/2/Pi(j).vnorm(k);
                Pd.ion(j).fv(k, :) = Pd.ion(j).fv(k, :)*numelPi(j)*Pi(j).vnorm(k);
                Pd.ion(j).fv0(k, :) = exp(-Pd.ion(j).v(k, :).^2)/sqrt(pi);
				Pd.ion(j).fv(k, :) = Pd.ion(j).fv(k, :) + Pd.ion(j).fv0(k, :);
            else
                [Pd.ion(j).fv(k, :), vsam(k, :)] = histcounts(P.ion(j).v(k, :), Ndist, 'Normalization', 'pdf');
				Pd.ion(j).v(k, :) = (vsam(k, 1:end-1) + vsam(k, 2:end))/2;
                Pd.ion(j).fv0(k, :) = exp(-Pd.ion(j).v(k, :).^2)/sqrt(pi);
		        if nargin >= 3
		            Pd.ion(j).fv(k, :) = Pd.ion(j).fv(k, :)*Pi(j).vnorm(k);
		            Pd.ion(j).v(k, :) = Pd.ion(j).v(k, :)/Pi(j).vnorm(k);
		        end
            end
        end
        if ifdf
            [Pd.ion(j).fvperp, ~, ~] = histcountsdf2(P.ion(j).v(2, :), P.ion(j).v(1, :), P.ion(j).w, [Ndist, Ndist], vsam(2, :), vsam(1, :));
            Pd.ion(j).fvperp = Pd.ion(j).fvperp*Pi(j).vnorm(1)*Pi(j).vnorm(2)*numelPi(j);
            Pd.ion(j).fvperp0 = exp(-(Pd.ion(j).v(2, :).').^2 - Pd.ion(j).v(1, :).^2)/pi;
            Pd.ion(j).fvperp = Pd.ion(j).fvperp + Pd.ion(j).fvperp0;
        else
            [Pd.ion(j).fvperp, ~, ~] = histcounts2(P.ion(j).v(2, :), P.ion(j).v(1, :), vsam(2, :), vsam(1, :), 'Normalization', 'pdf');
            Pd.ion(j).fvperp = Pd.ion(j).fvperp*Pi(j).vnorm(1)*Pi(j).vnorm(2);
            Pd.ion(j).fvperp0 = exp(-(Pd.ion(j).v(2, :).').^2 - Pd.ion(j).v(1, :).^2)/pi;
        end
    end
end
