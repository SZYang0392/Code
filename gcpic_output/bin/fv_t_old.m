function varargout = fv_t_old(fnum, Ps, Pause, type, Ndist, tscale, vphi)
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
    %!!!!!!!!!!!!!!!!!!!!!Can hold only one species of electrons and ions.

    % ------------ Initialize ------------ %
    e = 1.6021766208e-19;
    % Default parameters
    if nargin < 7
        vphi = 0;
    end
    if nargin < 6
        tscale = -1;
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
    
    if isscalar(fnum) && fnum < 0
        fnum = 1:(-fnum);
    end

    NP = numel(Ps);
    if strcmp(type, 'log')
        iflog = 1;
    end

    for k = 1:NP
        Ps(k).vt = sqrt(2*Ps(k).T*e/Ps(k).m);
        Ps(k).res = vphi/Ps(k).vt;
    end

    % ------------ Initial distribution ------------ %
    fv_init = distp_init(Ps, 'Ndist', Ndist);
    if iflog
        for k = 1:NP
            fv_init(k).logfv = log(fv_init(k).fv);
        end
    end

    % ------------ Figure ------------ %
    subnote = {'x', 'y', 'z'};
    Fig = cell(NP, 3);
    Inifig = cell(NP, 3);
    for k = 1:NP
        for j = 1:3
            Fig{k, j} = figure;
            hold on;
            if iflog
                Inifig{k, j} = plot(fv_init(k).v(j, :), fv_init(k).logfv(j, :), 'k-', 'LineWidth', 1, 'DisplayName', 'Initial');
                ylabel(['ln $f(v_', subnote{j}, ')$'], 'Interpreter', 'Latex');
            else
                Inifig{k, j} = plot(fv_init(k).v(j, :), fv_init(k).fv(j, :), 'k-', 'LineWidth', 1, 'DisplayName', 'Initial');
                ylabel(['Normalized $f(v_', subnote{j}, ')$'], 'Interpreter', 'Latex');
            end
            xlim([-5, 5]);
            xline(Ps(k).res, 'b-', 'LineWidth', 1, 'DisplayName', 'Resonant +');
            xline(-Ps(k).res, 'b-', 'LineWidth', 1, 'DisplayName', 'Resonant -');
            xlabel(['$v_', subnote{j}, ' / v_t$'], 'Interpreter', 'Latex');
            legend('Interpreter', 'Latex', 'Location', 'South');
            title([Ps(k).Name, ' $v_', subnote{j}, '$ distribution'], 'Interpreter', 'Latex');
            set(gca, 'FontSize', 12);
        end
    end

    % ------------ Time evolution ------------ %
    % Initialize data
    data = cell(size(Ps));         % fields : v, fv; [time, v, dim];
    for k = 1:NP
        data{k}.v = zeros(numel(fnum), Ndist, 3);
        data{k}.fv = zeros(numel(fnum), Ndist, 3);
        if iflog
            data{k}.logfv = zeros(numel(fnum), Ndist, 3);
        end
    end

    % Main Data Process
    fprintf('Loading and diagnosing particles...\n');
    t = tscale*fnum;          % wci*t
    for l = 1:numel(fnum)
        % Load particle
        pfname_mat = ['particle', num2str(fnum(l), '%04d'), '.mat'];
        if exist(pfname_mat, 'file')
            load(pfname_mat);
            Pdata = P;
            clear('P');
        else
            Pdata = loadp_f(fnum(l));
        end
        % Check T
        if tscale == -1
            t(l) = Pdata.param.stime;
        end
        % Pdata = loadparticle(fnum(l));
        Pd = distp(Pdata, Ndist, fv_init);
        % Collect distribution
        data{1}.v(l, :, :) = reshape(Pd.ele.v.', 1, Ndist, 3);
        data{1}.fv(l, :, :) = reshape(Pd.ele.fv.', 1, Ndist, 3);
        data{2}.v(l, :, :) = reshape(Pd.ion.v.', 1, Ndist, 3);
        data{2}.fv(l, :, :) = reshape(Pd.ion.fv.', 1, Ndist, 3);
        if iflog
            data{1}.logfv = log(data{1}.fv);
            data{2}.logfv = log(data{2}.fv);
        end
        % Display Process
        fprintf('%d / %d\n', l, numel(fnum));
    end

    % Figure folder
    dirname = 'Vdist';
    if ~exist(dirname, 'dir');
        mkdir(dirname);
    end

    % Main Figure Process
    Evofig = cell(NP, 3);
    if Pause
        I = cell(NP, 3);
        map = cell(NP, 3);
        for l = 1:numel(fnum)
            Lgd = ['$\omega_{ci}t = ', num2str(t(l), '%07f'), '$'];
            for k = 1:NP
                for j = 1:3
                    figure(Fig{k, j});
                    filename = [dirname, '\', Ps(k).Name, 'v', subnote{j}, '.gif'];
                    if l == 1
                        if iflog
                            Evofig{k, j} = plot(data{k}.v(l, :, j), data{k}.logfv(l, :, j), 'LineWidth', 1, 'DisplayName', Lgd);
                        else
                            Evofig{k, j} = plot(data{k}.v(l, :, j), data{k}.fv(l, :, j), 'LineWidth', 1, 'DisplayName', Lgd);
                        end
                        F = getframe(gcf);
                        im = frame2im(F);
                        [I{k, j}, map{k, j}] = rgb2ind(im,256);
                        imwrite(I{k, j}, map{k, j}, filename, 'GIF', 'Loopcount', inf, 'DelayTime', Pause);
                    else
                        if iflog
                            set(Evofig{k, j}, 'xdata', data{k}.v(l, :, j), 'ydata', data{k}.logfv(l, :, j), 'DisplayName', Lgd);
                        else
                            set(Evofig{k, j}, 'xdata', data{k}.v(l, :, j), 'ydata', data{k}.fv(l, :, j), 'DisplayName', Lgd);
                        end
                        F = getframe(gcf);
                        im = frame2im(F);
                        [I{k, j}, map{k, j}] = rgb2ind(im,256);
                        imwrite(I{k, j}, map{k, j}, filename, 'GIF', 'WriteMode', 'append', 'DelayTime', Pause);
                    end
                end
            end
        end
    else
        for l = 1:numel(fnum)
            Lgd = ['$\omega_{ci}t = ', num2str(t(l), '%07f'), '$'];
            for k = 1:NP
                for j = 1:3
                    figure(Fig{k, j});
                    if iflog
                        plot(data{k}.v(l, :, j), data{k}.logfv(l, :, j), 'LineWidth', 1, 'DisplayName', Lgd);
                    else
                        plot(data{k}.v(l, :, j), data{k}.fv(l, :, j), 'LineWidth', 1, 'DisplayName', Lgd);
                    end
                end
            end
        end
        % Save Figure;
        for k = 1:NP
            for j = 1:3
                filename = [dirname, '\', Ps(k).Name, 'v', subnote{j}];
                savefig(Fig{k, j}, [filename, '.fig']);
                saveas(Fig{k, j}, [filename, '.png']);
            end
        end
    end

    % Output
    if nargout == 1
        varargout{1} = data;
        for k = 1:NP
            for j = 1:3
                close(F{k, j});
            end
        end
    elseif nargout > 1
        varargout{2} = Fig;
    end
end