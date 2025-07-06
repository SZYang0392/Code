function varargout = XVspace_t(fnum, Ps, Pause, type, Ndist, tscale)
    % Show temporal evolution of particle distribution function.
    % ------------ Output ------------ %
    % varargout = [data, Fig];
    % Fig : Output figure
    % data : struct with fields : v - sample points (normalized to thermal velocity).
    %                             fv - distribution function.
    % ------------ Input ------------ %
    % fnum : filerange. Scalar: 1:fnum                                                 ----- compulsory
    %                   Vector: fnum(1), fnum(2), ... , fnum(end)
    % Ps : Particle struct VECTOR with at least fields 'vt'(or 'T' & 'm') & 'Name'.    ----- compulsory
    % type : 'log' if display ln f(v), otherwise display f(v).
    % Ndist : number of sample points
    % Pause : time(s) of changing figure. 0 if not change.
    % tscale : fnum*tscale = wci*t; tscale = dt*me/mi*ndiagp

    % ！！！!diagnosis_yang.f90里面按照fraci来区分粒子种类，换言之fraci相同的粒子种类无法区分

    % ------------ Initialize ------------ %
    % Default parameters
    if nargin < 6
        tscale = 1;
    end
    if nargin < 5
        Ndist = 150;
    end
    if nargin < 4
        type = 'log';
    end
    if nargin < 3
        Pause = 0.3;
    end
    
    if isscalar(fnum)
        fnum = 1:fnum;
    end

    NP = numel(Ps);
    for k = 1:NP
        if ~isfield(Ps(k), 'vt')
            Ps(k).vt = sqrt(2*Ps(k).T*e/Ps(k).m);
        elseif isempty(Ps(k).vt)
            Ps(k).vt = sqrt(2*Ps(k).T*e/Ps(k).m);
        end
    end

    % ------------ Initialize ------------ %
end