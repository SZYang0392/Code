function imagF_mp4(Find, rname, field, FFT, Unit, Frate)
    % Show temporal evolution of particle distribution function.
    % ------------ Input ------------ %
    % Find : index of file to be analyzed
    % rname : 2-element cell array, the direction to be plot, 'x', 'y', or 'z'
    % field : field of struct created in "fftkfield.m"
    % Unit : the unit of fourier spectrum, in LaTeX

    % ------------ Initialize ------------ %
    % Default parameters
    if nargin < 6
        Frate = 5;
    end
    if nargin < 5
        Unit = 'a.u.';
    end
    label = iscell(Unit) && numel(Unit) == 3;
    if nargin < 4
        FFT = false;
    end

    % ------------ Create Video ------------ %
    VideoName = [field, '2D.avi'];
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
        if FFT
            F = fftkfield(F);
        end
        Xdata = F.(rname{1});
        Ydata = F.(rname{2});
        Zdata = F.(field);
        if ~isreal(Zdata)
            Zdata = abs(Zdata);
        end

        % ------------ Plot ------------ %
        figure(Ffig);
        if k == 1
            Fimag = imagesc(Xdata, Ydata, Zdata);
            colorbar;
            set(gca, 'YDir', 'Normal');
            set(gca, 'Fontsize', 13);
            if label
                xlabel(['$', rname{1}, '\quad/\quad ', Unit{1}, '$'], 'Interpreter', 'latex');
                ylabel(['$', rname{2}, '\quad/\quad ', Unit{2}, '$'], 'Interpreter', 'latex');
            else
                if FFT
                    xlabel(['$k_', rname{1}, '\quad/\quad m^{-1}$'], 'Interpreter', 'Latex');
                    ylabel(['$k_', rname{2}, '\quad/\quad m^{-1}$'], 'Interpreter', 'Latex');
                else
                    xlabel(['$', rname{1}, '\quad/\quad m$'], 'Interpreter', 'Latex');
                    ylabel(['$', rname{2}, '\quad/\quad m$'], 'Interpreter', 'Latex');
                end
            end
        else
            set(Fimag, 'CData', Zdata);
        end

        if label
            title(['$', field, '\quad/\quad ', Unit{3}, '$ ', ' at $\omega_{ci}t=', num2str(F.stime), '$'], 'Interpreter', 'Latex');
        else
            title(['$', field, '\quad/\quad ', Unit, '$ ', ' at $\omega_{ci}t=', num2str(F.stime), '$'], 'Interpreter', 'Latex');
        end
        frame = getframe(gcf);
        writeVideo(VF, frame);
    end

    close(VF);
end