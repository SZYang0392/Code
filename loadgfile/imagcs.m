function F = imagcs(R, Z, frz, Opts)
    arguments
        R
        Z
        frz
        Opts.F
        Opts.Clear0 = false
        Opts.Ncolor = 20
        Opts.Alpha = 0.7
        Opts.FontSize = 13
        Opts.Clabel = 'a.u.'
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

    % Set alpha data - 1
    alphaOK = ~isnan(frz) & ~isinf(frz);
    if Opts.Clear0
        alphaOK = alphaOK & (frz~=0);
    end
    % Plot
    H = imagesc(R, Z, frz, 'AlphaData', alphaOK);
    set(gca, 'YDir', 'normal');
    axis equal;
    % Set alpha data - 1
    H.AlphaData = repmat(Opts.Alpha, size(frz));
    H.AlphaData(~alphaOK) = 0;
    % Set color
    colormap(jet(Opts.Ncolor));
    ax = gca;
    Color = clim(ax);
    clim([Color(1) - 0.6*(Color(2)-Color(1)), Color(2)]);
    % Set colorbar
    C = colorbar;
    C.Limits = [Color(1), C.Limits(2)];
    C.Label.String = Opts.Clabel;
    C.Label.Interpreter = 'latex';
    % Set legend and fontsize
    hold on;
    xlabel('$R/m$', 'Interpreter', 'latex');
    ylabel('$Z/m$', 'Interpreter', 'latex');
    set(gca, 'FontSize', Opts.FontSize);
end