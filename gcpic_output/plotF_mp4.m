function plotF_mp4(Find, rname, field, Unit, Frate)
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
    Min0 = 0;
    for k = 1:numel(Find)
        F = loadfield(Find(k));
        if ~F.OK
            continue;
        end
        if k == 1
            Max0 = max(F.(field), [], 'all');
            Max0 = max(Max0, 0);
            Min0 = min(F.(field), [], 'all');
            Min0 = min(Min0, 0);
        else
            Max1 = max(F.(field), [], 'all');
            Max0 = max(Max0, Max1);
            Min1 = min(F.(field), [], 'all');
            Min0 = min(Min0, Min1);
        end
    end


    % ------------ Create Video ------------ %
    VideoName = [field, '.avi'];
    VF = VideoWriter(VideoName, 'Motion JPEG AVI');
    VF.FrameRate = Frate;
    % VF.Quality = 95;
    open(VF);
    Ffig = figure;

    for k = 1:numel(Find)
        % ------------ Load and process data ------------ %
        F = loadfield(Find(k));
        if ~F.OK
            continue;
        end
        Xdata = F.(rname);
        Ydata = F.(field);

        % ------------ Plot ------------ %
        figure(Ffig);
        if k == 1
            Fcurve = plot(Xdata, Ydata, 'LineWidth', 1);
            set(gca, 'YDir', 'Normal');
            set(gca, 'Fontsize', 13);
            if Min0 <= 0 && Max0 >= 0
                ylim([1.1*Min0, 1.1*Max0]);
            elseif Min0 > 0 && Max0 >= 0
                ylim([0.9*Min0, 1.1*Max0]);
            elseif Min0 <= 0 && Max0 < 0
                ylim([1.1*Min0, 0.9*Max0]);
            end
            % ylim([-Inf, Inf]);
            xlabel(['$', rname, '\quad/\quad m^{-1}$'], 'Interpreter', 'Latex');
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