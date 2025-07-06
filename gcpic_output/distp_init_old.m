function Pd = distp_init_old(Ps, Ndist, vnorm, vmax)
% Return the initial Maxwellian distribution of velocity.
    % Ps is a particle struct array with fields : m, T
    % Ndist : number of sample points.
    % vnorm : normalization factor (vt, in m/s), basically nonzero
    %         vnorm = (0, 1, 0) denotes that normalization factors are set as sqrt(T_perp, T_para, T_perp)/sqrt(m).

    % --------------- Some default input parameters --------------- %
    if nargin < 2
        Ndist = 80;
    end
    if nargin < 4
        vmax = 3.5;                 % in normalization factor
    end

    % --------------- Pre calculation --------------- %
    e = 1.6021766208e-19;
    Pd = Ps;
    for k = 1:numel(Pd)
        if ~isfield(Pd(k), 'vt')
            Pd(k).vt = sqrt(2*Pd(k).T*e/Pd(k).m);
        elseif isempty(Pd(k).vt)
            Pd(k).vt = sqrt(2*Pd(k).T*e/Pd(k).m);
        end
    end

    % --------------- Normalization factor --------------- %
    % Case 1 : Default normalization factor
    if nargin < 3
        for k = 1:numel(Pd)
            Pd(k).vnorm = repmat(Pd(k).vt, 1, 3);
        end
    end
    if nargin >= 3
        % Case 2 : Universal normalization factor
        if isscalar(vnorm)
            for k = 1:numel(Pd)
                Pd(k).vnorm = repmat(vnorm, 1, 3);
            end
        end
        % Case 3 : Normalized to parallel and perpendicular thermal velocities
        if all(sort(vnorm) - [0, 0, 1] == 0, 'all')
            indpara = find(vnorm == 1);
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
            end
        end
    end

    % --------------- Normalized distribution function --------------- %
    sam = linspace(-vmax, vmax, Ndist);
    expsam = exp(-sam.^2);
    
    for k = 1:numel(Ps)
        Pd(k).v = zeros(3, Ndist);
        Pd(k).fv = zeros(3, Ndist);
        for j = 1:3
            Pd(k).v(j, :) = linspace(-vmax, vmax, Ndist);
            Pd(k).fv(j, :) = expsam/sqrt(pi);
        end
    end
end