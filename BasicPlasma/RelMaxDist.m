function varargout = RelMaxDist(m, T, Opts)
% Another set of input : Ps, v
% Calculate the normalized relativistic distribution.
% The momentum p is normalized to mc, vis normalized to c..
arguments
    m
    T = []
    Opts.p
    Opts.pmax
    Opts.pmN = 3.5
    Opts.dpN = 4
end
    % Physical constants
    c = 299792458;
    e = 1.602176565e-19;
    % Check input type
    if isstruct(m)
        Ps = m;
        m = Ps.m;
        T = Ps.T;
    end
    T = T(:).';
    % Determine the momentum array
    E02T = m*c^2/e./T;
    E02TMax = max(E02T, [], 'all');
    E02TMin = min(E02T, [], 'all');
    if isfield(Opts, 'p')
        p = Opts.p;
    else
        if isfield(Opts, 'pmax')
            pmax = Opts.pmax;
        else
            % Modity pmN here to adjust the maximum momentum
            pmax = sqrt((1 + Opts.pmN^2/E02TMin)^2 - 1);
        end
        % Modify dpN to adjust the spacing of the momentum array
        dp = sqrt((1 + 1./E02TMax).^2 - 1)/Opts.dpN;
        Np = ceil(pmax/dp) + 1;
        p = linspace(0, pmax, Np).';
    end

    % Calculate the normalization factor
    fNorm = E02T/(4*pi)./besselk(2, E02T, 1)./exp(-E02T);
    fp3d = fNorm.*exp(-E02T.*sqrt(1 + p.^2));
    fpmod = fp3d*(4*pi).*p.^2;
    v = 1./sqrt(1 + p.^(-2));
    E = sqrt(1 + p.^2) - 1;

    % The output
    varargout{1} = p;
    varargout{2} = fp3d;
    varargout{3} = fpmod;
    varargout{4} = v;
    varargout{5} = E;
end