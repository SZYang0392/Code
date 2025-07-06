function F = plotray(gen, G, rayf, Opts)
    arguments
        gen
        G = false
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
        Opts.Ind = 1:size(gen.(rayf), 2);
    end

    % Plot rays
    Rayf = gen.(rayf)(:, Opts.Ind, Opts.layer);
    if Opts.Normalize
        Rayf = Rayf./Rayf(1, :);
    end
    hold on;
    patch(gen.wr(:, Opts.Ind), gen.wz(:, Opts.Ind), Rayf, 'DisplayName', 'Rays', ...
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
    if ~isa(G, 'logical')
        hold on;
        plot(G.Raxis, G.Zaxis, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k', 'DisplayName', 'Axis');
        plot(G.Rbound, G.Zbound, 'r-', 'LineWidth', 1, 'DisplayName', 'LCFS');
        plot(G.Rlimiter, G.Zlimiter, 'k-', 'LineWidth', 1, 'DisplayName', 'Limiter');
        if Opts.Setlim
            Rmin = min(G.Rlimiter);
            Rmax = max(G.Rlimiter);
            dR = 0.15*(Rmax - Rmin);
            xlim([Rmin - dR, Rmax + dR]);
            Zmin = min(G.Zlimiter);
            Zmax = max(G.Zlimiter);
            dZ = 0.1*(Zmax - Zmin);
            ylim([Zmin - dZ, Zmax + dZ]);
        end
    end

    % Plot initial points
    if Opts.Plotinit
        hold on;
        plot(gen.wr(1, Opts.Ind), gen.wz(1, Opts.Ind), 'bo', 'MarkerSize', 2, 'MarkerFaceColor', 'b', ....
            'DisplayName', 'Launch');
    end

    set(gca, 'FontSize', Opts.FontSize);
end