function F = plotrayphase(genrayout, rayf, Opts)
    arguments
        genrayout
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
    if all(genrayout.wnpar(1, :) <= 0, 'all')
        genrayout.wnpar = -genrayout.wnpar;
        Ylabel = '$-N_{\parallel}$';
    else
        Ylabel = '$N_{\parallel}$';
    end
    patch(genrayout.spsi(:, Opts.Ind), genrayout.wnpar(:, Opts.Ind), Rayf, 'DisplayName', 'Rays', ...
    'edgecolor', 'flat', 'facecolor', 'none', 'MarkerFaceColor', 'flat');
    % Set figure
    xlabel('$\rho$', 'Interpreter', 'latex');
    ylabel(Ylabel, 'Interpreter', 'latex');
    title(Opts.title, 'Interpreter', 'latex');
    % Set color
    colormap(jet(Opts.Ncolor));
    C = colorbar;
    C.Label.String = Opts.unit;
    C.Label.Interpreter = 'latex';
    % Set xlim
    if Opts.Setlim
        Xmin = min(genrayout.spsi(:, Opts.Ind), [], 'all');
        Xmax = max(genrayout.spsi(:, Opts.Ind), [], 'all');
        dX = (Xmax - Xmin)*0.1;
        xlim([Xmin - 2*dX, Xmax + 0.2*dX]);
    end

    % Plot initial points
    if Opts.Plotinit
        hold on;
        plot(genrayout.spsi(1, Opts.Ind), genrayout.wnpar(1, Opts.Ind), 'co', 'MarkerSize', 5, ...
            'MarkerFaceColor', 'c', 'DisplayName','Launch');
    end

    set(gca, 'FontSize', Opts.FontSize);
end