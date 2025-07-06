function varargout = XVLeiQi_mp4(fnum, Ps, dim, xrange, mx, xshift, Ndist, vphi)
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
    % vpi : phase velocity of waves / wci

    % ！！！!diagnosis_yang.f90里面按照fraci来区分粒子种类，换言之fraci相同的粒子种类无法区分

    e = 1.6021766208e-19;

    % ------------ Initialize ------------ %
    % Default parameters
    if nargin < 8
        vphi = 0;
    end
    if nargin < 7
        Ndist = [500, 500];
    end
    if nargin < 6
        xshift = 0;
    end

    NP = numel(Ps);
    for k = 1:NP
        Ps(k).vt = sqrt(2*Ps(k).T*e/Ps(k).m);
    end

    % ------------ Create Video ------------ %
    Subscript = {'x', 'y', 'z'};
    videonameE = ['Ele', Subscript{dim}, '_', 'V', Subscript{dim}, '.avi'];
    videonameI = ['Ion', Subscript{dim}, '_', 'V', Subscript{dim}, '.avi'];
    VE = VideoWriter(videonameE, 'Motion JPEG AVI');
    VE.FrameRate = 5;
    % VE.Quality = 95;
    open(VE);
    VI = VideoWriter(videonameI, 'Motion JPEG AVI');
    VI.FrameRate = 5;
    % VI.Quality = 95;
    open(VI);
    FE = figure;
    FI = figure;

    for k = 1:fnum
        % ------------ Load and process data ------------ %
        filename = ['particle', num2str(k, '%04d'), '.mat'];
        if exist(filename, 'file')
            load(filename);
        else
            P = loadp_f(k);
        end
        if ~P.OK
            continue;
        end

        Ndistx = Ndist(1);
        Ndistv = Ndist(2);
        vmesh = linspace(1.5, 5.5, Ndistv + 1);
        xmesh = linspace(0, xrange/mx, Ndistx + 1);
        
        % P.ele.v(dim, :) = P.ele.v(dim, :)/Ps(1).vt;
        % % P.ele.r(dim, :) = mod(P.ele.r(dim, :) - xshift, (xrange/mx));
        % P.ion.v(dim, :) = P.ion.v(dim, :)/Ps(2).vt;
        % % P.ion.r(dim, :) = mod(P.ion.r(dim, :) - xshift, (xrange/mx));


        xshift1 = xshift + vphi*P.param.stime;
        P.ele.v([2, 3], :) = P.ele.v([2, 3], :)/Ps(1).vt;
        P.ele.r(dim, :) = mod(P.ele.r(dim, :) - xshift1, (xrange/mx));
        P.ion.v([2, 3], :) = P.ion.v([2, 3], :)/Ps(2).vt;
        P.ion.r(dim, :) = mod(P.ion.r(dim, :) - xshift1, (xrange/mx));

        [Ne0, Ve, Xe] = histcounts2(sqrt(sum(P.ele.v([2, 3], :).^2, 1)), P.ele.r(dim, :), vmesh, xmesh, 'Normalization', 'pdf');
        Xe = (Xe(2:end) + Xe(1:end-1))/2;
        Ve = (Ve(2:end) + Ve(1:end-1))/2;
        Ne = movmean(Ne0, 51, 1);
        Ne = movmean(Ne, 51, 2);
        [Ni0, Vi, Xi] = histcounts2(sqrt(sum(P.ion.v([2, 3], :).^2, 1)), P.ion.r(dim, :), vmesh, xmesh, 'Normalization', 'pdf');
        Xi = (Xi(2:end) + Xi(1:end-1))/2;
        Vi = (Vi(2:end) + Vi(1:end-1))/2;
        Ni = movmean(Ni0, 51, 1);
        Ni = movmean(Ni, 51, 2);

        % ------------ Plot ------------ %
        figure(FE);
        imagesc(Xe, Ve.', log(Ne));
        % ylim([1.5, 4]);
        colorbar;

        % plot(P.ele.r(dim, :), P.ele.v(dim, :), 'k.', 'Markersize', 0.000000005);
        % ylim([-4, 4]);
        % plot(P.ele.r(dim, :), sqrt(sum(P.ele.v([2, 3], :).^2, 1)), 'k.', 'Markersize', 0.000000005);
        % ylim([0, 4]);
        set(gca, 'YDir', 'Normal');
        set(gca, 'Fontsize', 13);
        xlabel(['$', Subscript{dim}, '/m$'], 'Interpreter', 'Latex');
        % ylabel(['$v_', Subscript{dim}, '/v_t$'], 'Interpreter', 'Latex');
        ylabel(['$v_\perp', '/v_t$'], 'Interpreter', 'Latex');
        title(['Electron phase space at $\omega_{ci}t=', num2str(P.param.stime, '%06f'), '$'], 'Interpreter', 'Latex');
        frame = getframe(gcf);
        writeVideo(VE, frame);

        figure(FI);
        imagesc(Xi, Vi.', log(Ni));
        % ylim([1.5, 4]);
        colorbar;

        % plot(P.ion.r(dim, :), P.ion.v(dim, :), 'k.', 'Markersize', 0.000000005);
        % ylim([-4, 4]);
        % plot(P.ion.r(dim, :), sqrt(sum(P.ion.v([2, 3], :).^2, 1)), 'k.', 'Markersize', 0.000000005);
        % ylim([0, 4]);
        set(gca, 'YDir', 'Normal');
        set(gca, 'Fontsize', 13);
        xlabel(['$', Subscript{dim}, '/m$'], 'Interpreter', 'Latex');
        % ylabel(['$v_', Subscript{dim}, '/v_t$'], 'Interpreter', 'Latex');
        ylabel(['$v_\perp', '/v_t$'], 'Interpreter', 'Latex');
        title(['Ion phase space at $\omega_{ci}t=', num2str(P.param.stime, '%06f'), '$'], 'Interpreter', 'Latex');
        frame = getframe(gcf);
        writeVideo(VI, frame);
    end

    close(VE);
    close(VI);
    

    Pd.ele.f = Ne;
    Pd.ele.x = Xe;
    Pd.ele.v = Ve;
    Pd.ion.f = Ni;
    Pd.ion.x = Xi;
    Pd.ion.v = Vi;

    varargout{1} = Pd;
    varargout{2} = FE;
    varargout{3} = FI;
end