function plotFk(Find, field, K_arr, Unit, Gamma)
    if nargin < 4
        Unit = 'a.u.';
    end
    FIELD = zeros(size(Find));
    t = zeros(size(Find));
    for k = 1:numel(Find)
        F = loadfield(Find(k));
        if ~F.OK
            fprintf(['field', num2str(Find(k), '%05d'), '.dat do not exist\n']);
            continue;
        end
        t(k) = F.stime;
        Fk = fftkfield(F);

        if nargin < 3
            if k == 1
                [FIELD(k), IND] = max(abs(Fk.(field)), [], [1, 2, 3], 'linear');
            else
                F0 = Fk.(field);
                FIELD(k) = abs(F0(IND));
            end
        else
            IND = [1, 1, 1];
            if F.nx > 1
                [kx, IND(1)] = min(abs(Fk.kx - K_arr(1)));
            end
            if F.ny > 1
                [ky, IND(2)] = min(abs(Fk.ky - K_arr(2)));
            end
            if F.nz > 1
                [kz, IND(3)] = min(abs(Fk.kz - K_arr(3)));
            end
            F0 = Fk.(field);
            FIELD(k) = abs(F0(IND(1), IND(2), IND(3)));
        end
    end
    Ffield = figure;
    plot(t, log(FIELD/FIELD(1)), 'LineWidth', 1);
    if nargin >= 5
        hold on;
        plot(t, -sqrt(2)/2 - Gamma*t, 'LineWidth',1);
    end
    xlabel('$\omega_{ci}t$', 'Interpreter', 'Latex');
    ylabel(['$', Unit, '$'], 'Interpreter', 'Latex');
    title(['Spacial fourier transform of', field], 'Interpreter', 'Latex');
    set(gca, 'Fontsize', 13);
    if ~exist('plot', 'dir')
        mkdir('plot');
    end
    savefig(Ffield, ['plot/', field, '.fig']);
    saveas(Ffield, ['plot/', field, '.png']);
    % close(Ffield);
end