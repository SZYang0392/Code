function F = plotray(genrayout, Gfile, rayf, Opts)
    arguments
        genrayout
        Gfile = false
        rayf = 'delpwr'
        Opts.F
        Opts.Ind
        Opts.layer = 1
        Opts.Normalize = false
        Opts.Ncolor = 50
        Opts.unit = false
        Opts.title = ''
        Opts.FontSize = 13
        Opts.Plotinit = true
        Opts.Setlim = true
    end

    % Open axis
    if isfield(Opts, 'F')
        F = Opts.F;
        if isa(F, 'matlab.ui.Figure')
            figure(F);
        elseif isa(F, 'matlab.graphics.axis.Axes')
            axes(F);
        end
    else
        F = figure;
    end

    % Select Unit
    if ~isa(Opts.unit, 'char')
        Opts.unit = GENRAY.rayunit(rayf);
    end

    % Select which ray to be plotted
    if ~isfield(Opts, 'Ind')
        Opts.Ind = 1:size(genrayout.(rayf), 2);
    end

    % Plot rays
    Rayf = genrayout.(rayf)(:, Opts.Ind, Opts.layer);
    if Opts.Normalize
        Rayf = Rayf./Rayf(1, :);
    end
    hold on;
    patch(genrayout.wr(:, Opts.Ind), genrayout.wz(:, Opts.Ind), Rayf, 'DisplayName', 'Rays', ...
    'edgecolor', 'flat', 'facecolor', 'none', 'MarkerFaceColor', 'flat');
    % Set figure
    axis equal;
    xlabel('$R/m$', 'Interpreter', 'latex');
    ylabel('$Z/m$', 'Interpreter', 'latex');
    title(Opts.title, 'Interpreter', 'latex');
    % Set color
    colormap(jet(Opts.Ncolor));
    C = colorbar;
    C.Label.String = Opts.unit;
    C.Label.Interpreter = 'latex';

    % Plot LCFS, magnetic axis and limiter
    if ~isa(Gfile, 'logical')
        hold on;
        plot(Gfile.Raxis, Gfile.Zaxis, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k', 'DisplayName', 'Axis');
        plot(Gfile.Rbound, Gfile.Zbound, 'r-', 'LineWidth', 1, 'DisplayName', 'LCFS');
        plot(Gfile.Rlimiter, Gfile.Zlimiter, 'k-', 'LineWidth', 1, 'DisplayName', 'Limiter');
        if Opts.Setlim
            Rmin = min(Gfile.Rlimiter);
            Rmax = max(Gfile.Rlimiter);
            dR = 0.15*(Rmax - Rmin);
            xlim([Rmin - dR, Rmax + dR]);
            Zmin = min(Gfile.Zlimiter);
            Zmax = max(Gfile.Zlimiter);
            dZ = 0.1*(Zmax - Zmin);
            ylim([Zmin - dZ, Zmax + dZ]);
        end
    end

    % Plot initial points
    if Opts.Plotinit
        hold on;
        plot(genrayout.wr(1, Opts.Ind), genrayout.wz(1, Opts.Ind), 'bo', 'MarkerSize', 2, 'MarkerFaceColor', 'b', ....
            'DisplayName', 'Launch');
    end

    set(gca, 'FontSize', Opts.FontSize);
end