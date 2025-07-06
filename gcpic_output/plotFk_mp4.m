function plotFk_mp4(Find, rname, field, Unit, Frate, Xlim)
    % Show temporal evolution of particle distribution function.
    % ------------ Input ------------ %
    % Find : index of file to be analyzed
    % rname : the direction to be plot, 'x', 'y', or 'z'
    % field : field of struct created in "fftkfield.m"
    % Unit : the unit of fourier spectrum, in LaTeX

    % ------------ Initialize ------------ %
    % Default parameters
    if nargin < 5
        Frate = 5;
    end
    if nargin < 4
        Unit = 'a.u.';
    end

    % ------------ Find Maximun ------------ %
    Max0 = 0;
    for k = 1:numel(Find)
        F = loadfield(Find(k));
        if ~F.OK
            continue;
        end
        Fk = fftkfield(F);
        if k == 1
            Max0 = max(abs(Fk.(field)), [], 'all');
        else
            Max1 = max(abs(Fk.(field)), [], 'all');
            Max0 = max(Max0, Max1);
        end
    end


    % ------------ Create Video ------------ %
    VideoName = [field, 'k', '.avi'];
    VF = VideoWriter(VideoName, 'Motion JPEG AVI');
    VF.FrameRate = Frate;
    % VF.Quality = 95;
    open(VF);
    Ffig = figure;

    rfield = ['k', rname];
    for k = 1:numel(Find)
        % ------------ Load and process data ------------ %
        F = loadfield(Find(k));
        if ~F.OK
            continue;
        end
        Fk = fftkfield(F);
        Xdata = Fk.(rfield);
        Ydata = abs(Fk.(field));

        % ------------ Plot ------------ %
        figure(Ffig);
        if k == 1
            Fcurve = plot(Xdata, Ydata, 'LineWidth', 1);
            set(gca, 'YDir', 'Normal');
            set(gca, 'Fontsize', 13);
            if nargin > 5
                xlim(Xlim);
            end
            ylim([0, 1.1*Max0]);
            xlabel(['$', rfield, '\quad/\quad m^{-1}$'], 'Interpreter', 'Latex');
            ylabel(['$', field, '\quad/\quad ', Unit, '$'], 'Interpreter', 'Latex');
        else
            set(Fcurve, 'xdata', Xdata, 'ydata', Ydata);
        end

        title(['$', field, '$\quad ', 'at $\omega_{ci}t=', num2str(F.stime), '$'], 'Interpreter', 'Latex');
        frame = getframe(gcf);
        writeVideo(VF, frame);
    end

    close(VF);
end