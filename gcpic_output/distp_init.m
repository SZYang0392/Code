function Pd = distp_init(Ps, Opts)
% Return the initial Maxwellian distribution of velocity.
    % Ps is a particle struct array with fields : m, T
    % Opts.Ndist : number of sample points.
    % Opts.vnorm : normalization factor (vt, in m/s), basically nonzero
    %         Opts.vnorm = (0, 1, 0) denotes that normalization factors are set as sqrt(T_perp, T_para, T_perp)/sqrt(m).

    % --------------- Some default input parameters --------------- %
    arguments
        Ps struct
        Opts.Ndist (1, 1) {mustBeInteger, mustBeGreaterThanOrEqual(Opts.Ndist, 2)} = 80
        Opts.vnorm (1, 3) {mustBeReal, mustBeNonNan}
        Opts.vmax (1, 1) {mustBeReal, mustBeNonNan, mustBePositive} = 3.5  % in normalization factor
    end

    % --------------- Pre calculation --------------- %
    e = 1.6021766208e-19;
    Pd = Ps;
    for k = 1:numel(Pd)
        if ~isfield(Pd(k), 'vt')
            Pd(k).vt = repmat(sqrt(2*Pd(k).T*e/Pd(k).m), 1, 3);
        elseif isempty(Pd(k).vt)
            Pd(k).vt = repmat(sqrt(2*Pd(k).T*e/Pd(k).m), 1, 3);
        end
    end

    % --------------- Normalization factor --------------- %
    % Case 1 : Default normalization factor
    if ~isfield(Opts, 'vnorm')
        for k = 1:numel(Pd)
            Pd(k).vnorm = Pd(k).vt;
			Pd(k).vnormp = sqrt(mean(Pd(k).vnorm([1, 2]).^2));
        end
    else
        if all(sort(Opts.vnorm) - [0, 0, 1] == 0, 'all')
            % Case 2 : Normalized to parallel and perpendicular thermal velocities
            indpara = find(Opts.vnorm == 1);
            for k = 1:numel(Pd)
                if ~isfield(Pd(k), 'T_para')
                    Pd(k).T_para = Pd(k).T;
                elseif isempty(Pd(k).T_para)
                    Pd(k).T_para = Pd(k).T;
                end
                if ~isfield(Pd(k), 'T_perp')
                    Pd(k).T_perp = Pd(k).T;
                elseif isempty(Pd(k).T_perp)
                    Pd(k).T_perp = Pd(k).T;
                end
                if ~isfield(Pd(k), 'vt_para')
                    Pd(k).vt_para = sqrt(2*Pd(k).T_para*e/Pd(k).m);
                elseif isempty(Pd(k).vt_para)
                    Pd(k).vt_para = sqrt(2*Pd(k).T_para*e/Pd(k).m);
                end
                if ~isfield(Pd(k), 'vt_perp')
                    Pd(k).vt_perp = sqrt(2*Pd(k).T_perp*e/Pd(k).m);
                elseif isempty(Pd(k).vt_perp)
                    Pd(k).vt_perp = sqrt(2*Pd(k).T_perp*e/Pd(k).m);
                end
                Pd(k).vnorm = repmat(Pd(k).vt_perp, 1, 3);
                Pd(k).vnorm(indpara) = Pd(k).vt_para;
				Pd(k).vnormp = Pd(k).vt_perp;
            end
        else
            % Case 3 : Universal normalization factor for all species
			vnormp = sqrt(mean(Opts.vnorm([1, 2]).^2));
            for k = 1:numel(Pd)
                Pd(k).vnorm = Opts.vnorm;
				Pd(k).vnormp = vnormp;
            end
        end
    end

    % --------------- Normalized distribution function --------------- %
    sam = linspace(-Opts.vmax, Opts.vmax, Opts.Ndist);
    expsam = exp(-sam.^2)/sqrt(pi);
    expsam_perp = exp(- sam.^2 - (sam.').^2)/pi;
    
    for k = 1:numel(Ps)
        Pd(k).v = zeros(3, Opts.Ndist);
        Pd(k).fv = zeros(3, Opts.Ndist);
        for j = 1:3
            Pd(k).v(j, :) = sam;
            Pd(k).fv(j, :) = expsam;
        end
        Pd(k).fvperp = expsam_perp;
    end
end
