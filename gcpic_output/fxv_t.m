function varargout = fxv_t(fnum, Ps, Opts, plotOpts)
    % Show temporal evolution of particle distribution function.
    % ------------ Output ------------ %
    % varargout = [data, Fig];
    % Fig : Output figure
    % data : struct with fields : v - sample points (normalized to thermal velocity).
    %                             fv - distribution function.
    % ------------ Input ------------ %
    % fnum : filerange.                                                                                                 ----- compulsory
    % Ps : Particle struct VECTOR with at least fields 'T', 'm' & 'Name'.    ----- compulsory
    %       Electrons are placed at the beginning of Ps.
    %       If Opts.E is given, Ps must contain fields 'q' & 'm'.
    % Opts.iflog : = true if display ln f(v), otherwise display f(v).
    % Opts.Ndist : number of sample points
    % Opts.tscale : wci*t*tscale = w*t  =>  tscale = w/wci

    %% ------------ Default input parameters ------------ %
    arguments
        fnum (1, :) {mustBeInteger} = -1
        % ----------- Physical Parameters ----------- %
        Ps struct = [];
        % Specify B directions
        Opts.B (1, 3) {mustBeReal, mustBeNonNan}
        % Specify resonant regions
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ONLY ONE REGION AVAILABLE NOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Opts.w (:, 1) {mustBeReal, mustBePositive}
        Opts.ktot (:, 3) {mustBeReal, mustBeNonNan}
        Opts.E (:, 3) double
        % Maximum velocity / vt
        Opts.vmax (1, 1) {mustBeReal, mustBeNonNan, mustBeNonzero} = 3.5
        % ----------- Figure Parameters ----------- %
        % Frame rate
        Opts.frate (1, 1) {mustBePositive, mustBeInteger}= 8
        % If save each frame
        Opts.ifsaveframe {mustBeNumericOrLogical} = false
        % Whether use log plot
        Opts.iflog (1, 1) {mustBeNumericOrLogical} = true
        % Number of sample points
        Opts.Ndist (1, 2) {mustBeInteger, mustBeGreaterThanOrEqual(Opts.Ndist, 2)} = 60
        % Scale and label of time
        Opts.tscale (1, 1) {mustBeReal} = -1
        Opts.ttext string = '\omega_{ci}t'
        % Figure range
        Opts.vlim (1, 2) {mustBeReal, mustBeNonNan} = [-4, 4]
        Opts.perplim (2, 2) {mustBeReal, mustBeNonNan} = [-3.3, 3.3; -3.3, 3.3]
        % Color of resonant regions
        Opts.Color (1, :) cell = {'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'}                                                % Color of resonant regions
        Opts.FaceAlpha (1, 1) = 0.2                                                                                             % Parameter of transparency
        % Phase space options
        Opts.Bnorm (1, 1) double {mustBeReal, mustBeNonNan, mustBePositive}                         % |B| in T
        Opts.mSource (1, 3) {mustBeInteger, mustBePositive} = [1, 1, 1]                                     % Simlation range / wavelength
        Opts.SimuRegion (1, 3) {mustBeReal, mustBeNonNan} = [0, 0, 0]                                    % Simulation region
        Opts.rshift (1, 3) double {mustBeReal, mustBeNonNan} = [0, 0, 0]                                  % Initial position of viewer
        Opts.ifstats (1, 1) {mustBeNumericOrLogical} = false                                                     % If calculate xv space distribution
		% If delta-f simulation is used
		Opts.ifdf (1, 1) {mustBeNumericOrLogical} = false
        % If distribution is given
        Opts.data cell = {}
        % If close figure
        Opts.ifclose (1, 1) {mustBeNumericOrLogical} = false
        % Guess that this property controls MATLAB function plot
        plotOpts.?matlab.graphics.chart.primitive.Line
    end

    %% ========================== Initialize ========================== %
    e = 1.6021766208e-19;
    u = 1.660538921e-27;
    % me = 9.10938291e-31;
    if fnum(1) == -1
        Fname = dir('particle*.nc');
        if isempty(Fname)
            Fname = dir('particle*.mat');
        end
        if isempty(Fname)
            Fname = dir('particle*.dat');
        end
        fnum = 1:numel(Fname);
    end
    % ------------------------ Calculalte direction array ------------------------ %
    if isfield(Opts, 'B')
            Bdir = Opts.B/sqrt(sum(Opts.B.^2));
            if all(Bdir == [1, 0, 0], 'all')
                Bdiry = [0, 1, 0];
            else
                Bdiry = cross(Bdir, [1, 0, 0]);
                Bdiry = Bdiry./sqrt(sum(Bdiry.^2));
            end
            Bdirx = cross(Bdiry, Bdir);
            Bdirx = Bdirx./sqrt(sum(Bdirx.^2));
            Dir = [Bdirx; Bdiry; Bdir].';
    else
        Dir = eye(3);
    end
    % ------------------------ Adjust direction array if ktot is given ------------------------ %
    ifw = isfield(Opts, 'w');
    ifktot = isfield(Opts, 'ktot');
    if ifw ~= ifktot
        error('w & ktot donot match!');
    end
    if ifktot
        if numel(Opts.w) ~= size(Opts.ktot, 1)
            error('w & ktot donot match!');
        end
        if isfield(Opts, 'B')
            % Adjust directions if the first k is not parallel to B_0
            kdir0 = Opts.ktot(1, :)./sum(Opts.ktot(1, :).^2);
            if ~all(Bdir == kdir0, 'all')
                Bdiry = cross(Bdir, kdir0);
                Bdiry = Bdiry/sqrt(sum(Bdiry.^2));
                Bdirx = cross(Bdiry, Bdir);
                Bdirx = Bdirx./sqrt(sum(Bdirx.^2));
                Dir = [Bdirx; Bdiry; Bdir].';
            end
        end
        ktot = Opts.ktot*Dir;      % [k_perp(1), B_0 x k_perp(1), B_0]
    end
    % ------------------------ Calculate velocity normalization factors ------------------------ %
    ifPs = ~isempty(Ps);
    if ~ifPs
        Ps = defaultPs(fnum(1), Dir);
    end
    NP = numel(Ps);
    for k = 1:NP
        if ifPs
            if isfield(Ps(k), 'uniso') && isfield(Opts, 'B')
                if Ps(k).uniso
                    Ps(k).vt = sqrt(2*[Ps(k).T_perp, Ps(k).T_perp, Ps(k).T_para]*e/Ps(k).m);
                else
                    Ps(k).vt = repmat(sqrt(2*Ps(k).T*e/Ps(k).m), 1, 3);
                end
            else
                Ps(k).vt = repmat(sqrt(2*Ps(k).T*e/Ps(k).m), 1, 3);
            end
        end
    end
    % ------------------------ If w & k is given, calculate resonant velocities ------------------------ %
    if ifw
        ktot20 = sum((Opts.ktot.').^2).';
        vphi0 = Opts.w.*Opts.ktot./ktot20;
        ktot2 = sum((ktot.').^2).';
        vphi = Opts.w./ktot;
        absk = sqrt(ktot2(1));
        for k = 1:NP
            Ps(k).res = vphi(:, [1, 1, 3])./Ps(k).vt;              % [perp, perp, para]
            Ps(k).res_para = ktot(:, end).*Ps(k).vt(end)./(ktot(:, 1).*Ps(k).vt(1));
        end
    end
    % ------------------------ If electric field is given, calculate rsonant half_width ------------------------ %
    if isfield(Opts, 'E')
        E = Opts.E*Dir;
        E1 = abs(E);
        E1([1, 2]) = sqrt(sum(E([1, 2]).^2)/2);
        ktot1 = ktot;
        ktot1([1, 2]) = sqrt(sum(ktot([1, 2]).^2)/2);
        for k = 1:NP
            Ps(k).resw = 2*sqrt(abs(Ps(k).q)*e*E1./Ps(k).m./ktot1)./Ps(k).vt;
        end
        kdir = ktot(1, :)./sum(ktot1(1, :).^2);
    end
    % ------------------------ Plot options ------------------------ %
    if ~isempty(fieldnames(plotOpts))
        plotCell = namedargs2cell(plotOpts);
    else
        plotCell.LineWidth = 1;
    end
    
    %% ========================== Initial distribution ========================== %%
    fv_init = distp_init(Ps, 'Ndist', Opts.Ndist(1), 'vmax', Opts.vmax);
    if Opts.iflog
        for k = 1:NP
            fv_init(k).logfv = log(fv_init(k).fv);
        end
    end
    
    %% ========================== Theoretical Figure at t = 0 ========================== %%
    if isfield(Opts, 'B')
        subnote = {'k\perp', 'B_{0}\times k', 'B_{0}'};
    else
        subnote = {'x', 'y', 'z'};
    end
    Fig = cell(NP, 1);
    Inifig = cell(NP, 4);
    Cbar = cell(2, NP);
    for k = 1:NP
        Fig{k} = figure;
        set(Fig{k}, 'OuterPosition', [1, 1, 1440, 900]);
        for j = 1:3
            % ------------------------ Figure 1~3 : v - f(v) Distribution in 3 directions ------------------------ %
            Fsub = subplot(2, 4, j);
            hold on;
            if Opts.iflog
                Inifig{k, j} = plot(fv_init(k).v(j, :), fv_init(k).logfv(j, :), 'k-', plotCell, 'DisplayName', 'Maxwellian');
                ylabel(['$\ln\left[v_{t', subnote{j}, '}\cdot f(v_{', subnote{j}, '})\right]$'], 'Interpreter', 'Latex');
            else
                Inifig{k, j} = plot(fv_init(k).v(j, :), fv_init(k).fv(j, :), 'k-', plotCell, 'DisplayName', 'Maxwellian');
                ylabel(['$v_{t', subnote{j}, '}\cdot f(v_{', subnote{j}, '})$'], 'Interpreter', 'Latex');
            end
            % Resonant region
            if ifw
                if isfield(Opts, 'E')
                    YDot = [Fsub.YLim(1), Fsub.YLim(2), Fsub.YLim(2), Fsub.YLim(1), Fsub.YLim(1)];
                end
                for l = 1:numel(Opts.w)
                    if isfield(Opts, 'E')
                        XDot = [-1, -1, 1, 1, -1]*Ps(k).resw(l, j);
                    end
                    if ~isnan(Ps(k).res(l, j)) && ~isinf(Ps(k).res(l, j))
                        if j < 3
                            xline(Ps(k).res(l, j), 'LineWidth', 1, 'DisplayName', ['Res ', num2str(l), '+']);
                            xline(-Ps(k).res(l, j), 'LineWidth', 1, 'DisplayName', ['Res ', num2str(l), '-']);
                            if isfield(Opts, 'E')
                                fill(XDot + Ps(k).res(l, j), YDot, Opts.Color{l}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', ['ResRange ', num2str(l), '+']);
                                fill(XDot - Ps(k).res(l, j), YDot, Opts.Color{l}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', ['ResRange ', num2str(l), '-']);
                            end
                        else
                            xline(Ps(k).res(l, j), 'LineWidth', 1, 'DisplayName', ['Res ', num2str(l)]);
                            if isfield(Opts, 'E')
                                fill(XDot + Ps(k).res(l, j), YDot, Opts.Color{l}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', ['ResRange ', num2str(l)]);
                            end
                        end
                    end
                end
            end
            xlim(Opts.vlim);
            xlabel(['$v_{', subnote{j}, '} / v_{t', subnote{j}, '}$'], 'Interpreter', 'Latex');
            legend('Interpreter', 'Latex', 'Location', 'South');
            title(join([' $v_{', subnote{j}, '}$ distribution']), 'Interpreter', 'Latex');
            set(gca, 'FontSize', 12);
            % ------------------------ Figure 5~7 : x - v - f(v) Distribution in 3 directions ------------------------ %
            subplot(2, 4, 4 + j);
            hold on;
            % Resonant region
            if ifw
                for l = 1:numel(Opts.w)
                    if ~isnan(Ps(k).res(l, j)) && ~isinf(Ps(k).res(l, j))
                        if j < 3
                            yline(Ps(k).res(l, j), 'LineWidth', 1, 'DisplayName', ['Res ', num2str(l), '+']);
                            yline(-Ps(k).res(l, j), 'LineWidth', 1, 'DisplayName', ['Res ', num2str(l), '-']);
                        else
                            yline(Ps(k).res(l, j), 'LineWidth', 1, 'DisplayName', ['Res ', num2str(l)]);
                        end
                    end
                end
            end
            ylim(Opts.vlim);
        end
        % ------------------------ Figure 4 : v_\perp - f(v_\perp) Distribution ------------------------ %
        Fsub = subplot(2, 4, 4);
        hold on;
        Inifig{k, 4} = imagesc(fv_init(k).v(1, :), fv_init(k).v(2, :), zeros(size(fv_init(k).fvperp)));
        if k == 1
            C = colormap;
            C = C(1, :);
        end
        F0 = Inifig{k, 4}.Parent;
        F0.Color = C;
        clear('F0');
        set(gca, 'YDir', 'Normal');
        Cbar{k} = colorbar('TickLabelInterpreter', 'latex', 'FontSize', 12);
        Cbar{k}.Label.Interpreter = 'latex';
        if Opts.iflog
            Cbar{k}.Label.String = '$\ln \left[v_{t\perp}^2 f(v_{\perp})\right] - \ln \left[v_{t\perp}^2 f_M(v_{\perp})\right]$';
        else
            Cbar{k}.Label.String = '$v_{t\perp}^2 \left[f(v_{\perp}) - f_M(v_{\perp})]$';
        end
        set(gca, 'FontSize', 12);
        if ifw
            if isfield(Opts, 'E')
                YDot = [Fsub.YLim(1), Fsub.YLim(2), Fsub.YLim(2), Fsub.YLim(1), Fsub.YLim(1)];
            end
            for l = 1:numel(Opts.w)
                if ~isnan(Ps(k).res(l, 1)) && ~isinf(Ps(k).res(l, 1))
                    xline(Ps(k).res(l, 1), 'LineWidth', 1, 'DisplayName', ['Res ', num2str(l)]);
                    if isfield(Opts, 'E')
                        XDot = [-1, -1, 1, 1, -1]*Ps(k).resw(l, 1);
                        fill(XDot + Ps(k).res(l, 1), YDot, Opts.Color{l}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', ['ResRange ', num2str(l)]);
                    end
                end
            end
        end
        % legend('Location', 'SouthWest', 'Interpreter', 'latex');
        xlim(Opts.perplim(1, :));
        ylim(Opts.perplim(2, :));
        xlabel(['$v_{t', subnote{1}, '}$'], 'Interpreter', 'Latex');
        ylabel(['$v_{t', subnote{2}, '}$'], 'Interpreter', 'Latex');
        title('$v_{\perp}$ Maxwellian', 'Interpreter', 'Latex');
        set(gca, 'FontSize', 12);
        % ------------------------ Figure 8 : k*r - v_\perp - f(v_\perp) Distribution ------------------------ %
        subplot(2, 4, 8);
        hold on;
        % Resonant region
        if ifw
            for l = 1:numel(Opts.w)
                if ~isnan(Ps(k).res(l, 1)) && ~isinf(Ps(k).res(l, 1))
                    yline(Ps(k).res(l, 1), 'LineWidth', 1, 'DisplayName', ['Res ', num2str(l)]);
                end
            end
        end
        ylim([0, Opts.vlim(2)]);
        sgtitle(Ps(k).Name, 'FontSize', 15);
    end
    
    %% ========================== Initialize Data Process ========================== %%
    if isfield(Opts, 'Bnorm')
        wci = Opts.Bnorm*e/u;
    else
        wci = Inf;
    end
    if ifw
        vphi0_wci = vphi0(1, :).'/wci;
    end
    Lam0 = (Opts.SimuRegion./Opts.mSource).';
    fprintf('Loading and diagnosing particles...\n');
    t = 1:numel(fnum);          % wci*t
    % ------------------------ Attribute memory for struct data ------------------------ %
    data = cell(size(Ps));         % fields : v, fv; [time, v, dim];
    for k = 1:NP
        data{k}.v = zeros(3, Opts.Ndist(1), numel(fnum));
        data{k}.fv = zeros(3, Opts.Ndist(1), numel(fnum));
        data{k}.fvperp = zeros(Opts.Ndist(1), Opts.Ndist(1), numel(fnum));
        if Opts.iflog
            data{k}.logfv = zeros(3, Opts.Ndist(1), numel(fnum));
        end
        data{k}.r = zeros(3, Opts.Ndist(2), numel(fnum));
        data{k}.fxv = zeros(Opts.Ndist(2), Opts.Ndist(1), 3, numel(fnum));
        data{k}.vperpnorm = 0;
        data{k}.vperp = zeros(numel(fnum), Opts.Ndist(1));
        data{k}.rperp = zeros(numel(fnum), Opts.Ndist(2));
        data{k}.fxv_perp = zeros(Opts.Ndist(1), Opts.Ndist(2), numel(fnum));
    end

    %% ========================== Main Loop for Data Collecting and Figure ========================== %
    % ------------------------ Figure cell array ------------------------ %
    Evofig = cell(NP, 3);
    Evofigxv = cell(NP, 3);
    Evofigperp = cell(NP, 1);
    % ------------------------ Generate videos or figures ------------------------ %
    dirname = 'Vdist';
    if ~exist(dirname, 'dir')
        mkdir(dirname);
    end
    if Opts.frate > 0
        VideoName = cell(NP, 1);
        VF = cell(NP, 1);
        for k = 1:NP
            if isfield(Opts, 'B')
                VideoName{k} = [Ps(k).Name, '_xv.avi'];
                VF{k} = VideoWriter([dirname, '/', VideoName{k}], 'Motion JPEG AVI');
            else
                VideoName{k} = [Ps(k).Name, '_xv_B.avi'];
                VF{k} = VideoWriter([dirname, '/', VideoName{k}], 'Motion JPEG AVI');
            end
            VF{k}.FrameRate = Opts.frate;
            open(VF{k});
        end
    end
    % ------------------------ Main loop ------------------------ %
    for l = 1:numel(fnum)
        tic;
        % ------------------------ Load particle ------------------------ %
        pfname_mat = ['particle', num2str(fnum(l), '%04d'), '.mat'];
        if exist(pfname_mat, 'file')
            load(pfname_mat, 'P');
            Pdata = P;
            clear('P');
        else
            Pdata = loadp_f(fnum(l));
        end
        % ------------------------ Check time ------------------------ %
        if Opts.tscale == -1
            t(l) = Pdata.param.stime;
        else
            t(l) = Pdata.stime*Opts.tscale;
        end
        % ------------------------ Change direction ------------------------ %
        for k = 1:numel(Pdata.ele)
            Pdata.ele(k).v = (Dir.')*Pdata.ele(k).v;
            if ifw
                Pdata.ele(k).r = mod(Pdata.ele(k).r - Opts.rshift.' - vphi0_wci*Pdata.param.stime, Lam0);
            end
            Pdata.ele(k).r = (Dir.')*Pdata.ele(k).r;
        end
        for k = 1:numel(Pdata.ion)
            Pdata.ion(k).v = (Dir.')*Pdata.ion(k).v;
            if ifw
                Pdata.ion(k).r = mod(Pdata.ion(k).r - Opts.rshift.' - vphi0_wci*Pdata.param.stime, Lam0);
            end
            Pdata.ion(k).r = (Dir.')*Pdata.ion(k).r;
        end
        Pdata1 = [Pdata.ele, Pdata.ion];
        % ------------------------ Calculate and collect velocity distribution ------------------------ %
        EN = numel(Pdata.ele);
        Pd = distp_v(Pdata, Opts.Ndist(1), fv_init, EN, Opts.ifdf);
        Pd1 = [Pd.ele, Pd.ion];
        for k = 1:NP
            data{k}.v(:, :, l) = Pd1(k).v;
            data{k}.fv(:, :, l) = Pd1(k).fv;
            if Opts.iflog
                data{k}.logfv(:, :, l) = log(data{k}.fv(:, :, l));
            end
            data{k}.fvperp(:, :, l) = Pd1(k).fvperp;
        end
        % ------------------------ Calculate and collect phase space distribution ------------------------ %
        if Opts.ifstats
            PD = distp_xv(Pdata, Opts.Ndist, fv_init, EN, Opts.ifdf, [1, 2], kdir);
            PD1 = [PD.ele, PD.ion];
            for k = 1:NP
                data{k}.r(:, :, l) = PD1(k).r;
                data{k}.fxv(:, :, :, l) = PD1(k).fxv;
                data{k}.vperpnorm = PD1(k).vperpnorm;
                data{k}.vperp(l, :) = PD1(k).vperp;
                data{k}.rperp(l, :) = PD1(k).rperp;
                data{k}.fxv_perp(:, :, l) = PD1(k).fxv_perp;
            end
        end
        % ------------------------ Display Process ------------------------ %
        fprintf('Particle : %d / %d\n', l, numel(fnum));

        Lgd = join(['$', Opts.ttext, ' = ', num2str(t(l), '%06f'), '$']);
        for k = 1:NP
            figure(Fig{k});
            for j = 1:3
                % ------------------------ Figure 1~3 : v - f(v) Distribution in 3 directions ------------------------ %
                subplot(2, 4, j);
                if l == 1
                    if Opts.iflog
                        Evofig{k, j} = plot(data{k}.v(j, :, l), data{k}.logfv(j, :, l), plotCell, 'DisplayName', Lgd);
                    else
                        Evofig{k, j} = plot(data{k}.v(j, :, l), data{k}.fv(j, :, l), plotCell, 'DisplayName', Lgd);
                    end
                else
                    if Opts.iflog
                        set(Evofig{k, j}, 'xdata', data{k}.v(j, :, l), 'ydata', data{k}.logfv(j, :, l), 'DisplayName', Lgd);
                    else
                        set(Evofig{k, j}, 'xdata', data{k}.v(j, :, l), 'ydata', data{k}.fv(j, :, l), 'DisplayName', Lgd);
                    end
                end
                % ------------------------ Figure 5~7 : x - v - f(v) Distribution in 3 directions ------------------------ %
                Fsub = subplot(2, 4, j + 4);
                if l == 1
                    if Opts.ifstats
                        Evofigxv{k, j} = imagesc(data{k}.r(j, :, l), data{k}.v(j, :, l).', data{k}.fxv(:, :, l));
                    else
                        Evofigxv{k, j} = plot(Pdata1(k).r(j, :), Pdata1(k).v(j, :)/Ps(k).vt(j), 'k.', 'Markersize', 0.000000005);
                    end
                    clear('F0');
                    xlabel(['$r_{', subnote{j}, '} / m$   ,   $m_{', subnote{j}, '} = ', num2str(Opts.mSource(j)), '$'], 'Interpreter', 'latex');
                    ylabel(['$v_{', subnote{j}, '} / v_{t', subnote{j}, '}$'], 'Interpreter', 'Latex');
                    title(join(['$r_{', subnote{j}, '} - v_{', subnote{j}, '}$ figure at ', Lgd]), 'Interpreter', 'latex');
                    % Resonant region
                    if ifw && isfield(Opts, 'E')
                        for m = 1:numel(Opts.w)
                            if isnan(Ps(k).res(m, j)) && isinf(Ps(k).res(m, j))
                                continue;
                            end
                            XDot = [Fsub.XLim(1), Fsub.XLim(1), Fsub.XLim(2), Fsub.XLim(2), Fsub.XLim(1)];
                            YDot = [-1, 1, 1, -1, -1]*Ps(k).resw(m, j);
                            if j < 3
                                fill(XDot, YDot + Ps(k).res(m, j), Opts.Color{m}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', ['ResRange ', num2str(m), '+']);
                                fill(XDot, YDot - Ps(k).res(m, j), Opts.Color{m}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', ['ResRange ', num2str(m), '-']);
                            else
                                fill(XDot, YDot + Ps(k).res(m, j), Opts.Color{m}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', ['ResRange ', num2str(m)]);
                            end
                        end
                    end
                    set(gca, 'FontSize', 12);
                else
                    if Opts.ifstats
                        Evofigxv{k, j}.XData = data{k}.r(j, :, l);
                        Evofigxv{k, j}.YData = data{k}.v(j, :, l);
                        Evofigxv{k, j}.CData = data{k}.fxv(:, :, l);
                    else
                        Evofigxv{k, j}.XData = Pdata1(k).r(j, :);
                        Evofigxv{k, j}.YData = Pdata1(k).v(j, :)/Ps(k).vt(j);
                    end
                end
                title(join(['$r_{', subnote{j}, '} - v_{', subnote{j}, '}$ figure at ', Lgd]), 'Interpreter', 'latex');
            end
            % ------------------------ Figure 4 : v_\perp - f(v_\perp) Distribution ------------------------ %
            subplot(2, 4, 4);
            if Opts.iflog
                Inifig{k, 4}.CData = log(data{k}.fvperp(:, :, l)) + log(pi) + data{k}.v(1, :, l).^2 + (data{k}.v(2, :, l).').^2;
            else
                Inifig{k, 4}.CData = data{k}.fvperp(:, :, l) - (1/pi)*exp(- data{k}.v(1, :, l).^2 - (data{k}.v(2, :, l).').^2);
            end
            Inifig{k, 4}.XData = data{k}.v(1, :, l);
            Inifig{k, 4}.YData = data{k}.v(2, :, l);
            title(join(['$v_{\perp}$ Distribution at $', Opts.ttext, ' = ', num2str(t(l), '%06f'), '$']), 'Interpreter', 'Latex');
            % ------------------------ Figure 8 : k*r - v_\perp - f(v_\perp) Distribution ------------------------ %
            Fsub = subplot(2, 4, 8);
            if l == 1
                if Opts.ifstats
                    Evofigperp{k} = imagesc(data{k}.rperp(l, :)*absk, data{k}.vperp(l, :).', data{k}.fxv_perp(:, :, l)/absk);
                    F0 = Evofigperp{k}.Parent;
                    F0.Color = C;
                    clear('F0');
                else
                    Evofigperp{k} = plot(ktot*Pdata1(k).r, sqrt(sum(Pdata1(k).v([1, 2], :).^2))/Ps(k).vt(1), 'k.', 'Markersize', 0.000000005);
                end
                xlabel('$k\cdot r$', 'Interpreter', 'latex');
                ylabel('$v_{\perp} / v_{t\perp}$', 'Interpreter', 'latex');
                title(join(['$v_{\perp}$ Distribution at', Lgd]), 'Interpreter', 'Latex');
                % Resonant region
                if ifw && isfield(Opts, 'E')
                    for m = 1:numel(Opts.w)
                        if isnan(Ps(k).res(m, 1)) && isinf(Ps(k).res(m, 1))
                            continue;
                        end
                        XDot = [Fsub.XLim(1), Fsub.XLim(1), Fsub.XLim(2), Fsub.XLim(2), Fsub.XLim(1)];
                        YDot = [-1, 1, 1, -1, -1]*Ps(k).resw(m, 1);
                        fill(XDot, YDot + Ps(k).res(m, 1), Opts.Color{m}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', ['ResRange ', num2str(m)]);
                    end
                end
                set(gca, 'FontSize', 12);
            else
                if Opts.ifstats
                    Evofigperp{k}.XData = data{k}.rperp(l, :)*absk;
                    Evofigperp{k}.YData = data{k}.vperp(l, :).';
                    Evofigperp{k}.CData = data{k}.fxv_perp(:, :, l)/absk;
                else
                    Evofigperp{k}.XData = ktot*Pdata1(k).r;
                    Evofigperp{k}.YData = sqrt(sum(Pdata1(k).v([1, 2], :).^2))/Ps(k).vt(1);
                end
            end
            title(join(['$v_{\perp}$ Distribution at', Lgd]), 'Interpreter', 'Latex');
            % ------------------------ Save Frame ------------------------ %
            if Opts.frate > 0
                F = getframe(Fig{k});
                writeVideo(VF{k}, F);
            end
            if Opts.ifsaveframe
                if isfield(Opts, 'B')
                    filename = [dirname, '/', Ps(k).Name, '_xv_B', num2str(fnum(l))];
                else
                    filename = [dirname, '/', Ps(k).Name, '_xv', num2str(fnum(l))];
                end
                savefig(Fig{k}, [filename, '.fig']);
                saveas(Fig{k}, [filename, '.png']);
            end
        end
        % ------------------------ Display Process ------------------------ %
        fprintf('Frame : %d / %d\n', l, numel(fnum));
        toc;
    end
    if Opts.frate > 0
        for k = 1:NP
            close(VF{k});
        end
    end

    %% ========================== Output ========================== %
    varargout{1} = data;
    varargout{2}.t = t;
    varargout{2}.tscale = Opts.tscale;
    varargout{2}.ttext = Opts.ttext;
    if Opts.ifclose
        for k = 1:NP
            close(Fig{k});
        end
        varargout{3} = {};
    else
        varargout{3} = Fig;
    end
end

%% ========================== Functions for automatically generating Ps array ========================== %
function Ps = defaultPs(fid, Dir)
	pfname_mat = ['particle', num2str(fid, '%04d'), '.mat'];
    if exist(pfname_mat, 'file')
        load(pfname_mat, 'P');
        Pdata = P;
        clear('P');
    else
        Pdata = loadp_f(fid);
    end
    Ps = struct('Name', [], 'vt', []);
    Ne = numel(Pdata.ele);
    Pe = repmat(Ps, 1, Ne);
    for j = 1:Ne
        Pe(j).Name = ['Ele', num2str(j)];
        Pe(j).vt = sqrt(2*mean((Pdata.ele(j).v.'*Dir).^2));
    end
    Ni = numel(Pdata.ion);
    Pi = repmat(Ps, 1, Ni);
    for j = 1:Ni
        Pi(j).Name = ['Ion', num2str(j)];
        Pi(j).vt = sqrt(2*mean((Pdata.ion(j).v.'*Dir).^2));
    end
    Ps = [Pe, Pi];
end
