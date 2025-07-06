function varargout = XVLeiQi(fnum, Ps, dim, xrange, mx, xshift, Ndist)
    % Show temporal evolution of particle distribution function.
    % ------------ Output ------------ %
    % varargout = [data, Fig];
    % Fig : Output figure
    % data : struct with fields : v - sample points (normalized to thermal velocity).
    %                             fv - distribution function.
    % ------------ Input ------------ %
    % fnum : index of file to be analyzed
    % Ps : Particle struct VECTOR with at least fields 'vt'(or 'T' & 'm') & 'Name'.    ----- compulsory
    % dim : the coordinate to be analyzed.
    % xrange : simulation domain
    % mx : number of waves contained in the simulation domain
    % xshift : positon of particles are left-shifted this range.
    % Ndist : number of sample points

    % ！！！!diagnosis_yang.f90里面按照fraci来区分粒子种类，换言之fraci相同的粒子种类无法区分

    e = 1.6021766208e-19;

    % ------------ Initialize ------------ %
    % Default parameters
    if nargin < 7
        Ndist = [100, 100];
    end
    if nargin < 6
        xshift = 0;
    end

    NP = numel(Ps);
    for k = 1:NP
        Ps(k).vt = sqrt(2*Ps(k).T*e/Ps(k).m);
    end

    % ------------ Load and process data ------------ %
    filename = ['particle', num2str(fnum, '%04d'), '.mat'];
    if exist(filename, 'file')
        load(filename);
    else
        P = loadp_f(fnum);
    end

    Ndistx = Ndist(1);
    Ndistv = Ndist(2);
    vmesh = linspace(1.5, 5.5, Ndistv + 1);
    xmesh = linspace(0, xrange/mx, Ndistx + 1);
    
    P.ele.v(dim, :) = P.ele.v(dim, :)/Ps(1).vt;
    % P.ele.r(dim, :) = mod(P.ele.r(dim, :) - xshift, (xrange/mx));
    P.ion.v(dim, :) = P.ion.v(dim, :)/Ps(2).vt;
    % P.ion.r(dim, :) = mod(P.ion.r(dim, :) - xshift, (xrange/mx));

    [Ne, Ve, Xe] = histcounts2(abs(P.ele.v(dim, :)), P.ele.r(dim, :), vmesh, xmesh, 'Normalization', 'pdf');
    Xe = (Xe(2:end) + Xe(1:end-1))/2;
    Ve = (Ve(2:end) + Ve(1:end-1))/2;
    [Ni, Vi, Xi] = histcounts2(abs(P.ion.v(dim, :)), P.ion.r(dim, :), vmesh, xmesh, 'Normalization', 'pdf');
    Xi = (Xi(2:end) + Xi(1:end-1))/2;
    Vi = (Vi(2:end) + Vi(1:end-1))/2;

    Pd.ele.f = Ne;
    Pd.ele.x = Xe;
    Pd.ele.v = Ve;
    Pd.ion.f = Ni;
    Pd.ion.x = Xi;
    Pd.ion.v = Vi;

    % ------------ Plot ------------ %
    Fe = figure;
    % imagesc(Xe, Ve.', log(Ne));
    plot(P.ele.r(dim, :), P.ele.v(dim, :), 'k.', 'Markersize', 0.000000005);
    set(gca, 'YDir', 'Normal');
    set(gca, 'Fontsize', 13);
    colorbar;
    xlabel('x/m');
    ylabel('v/v_t');
    title('Electron phase space');
    Fi = figure;
    % imagesc(Xi, Vi.', log(Ni));
    plot(P.ion.r(dim, :), P.ion.v(dim, :), 'k.', 'Markersize', 0.000000005);
    set(gca, 'YDir', 'Normal');
    set(gca, 'Fontsize', 13);
    colorbar;
    xlabel('x/m');
    ylabel('v/v_t');
    title('Ion phase space');
    
    varargout{1} = Pd;
    varargout{2} = Fe;
    varargout{3} = Fi;
end