function plotF(Find, field, dt, Unit)
% dt in wlh^{-1}
    for k = Find
        F = loadfield(k);
        if ~F.OK
            fprintf(['field', num2str(k, '%05d'), '.dat do not exist\n']);
            continue;
        end
        FIELD = F.(field);
        Ffield = figure;
        imagesc(F.x, F.y, FIELD);
        colorbar;
        xlabel('$x / m$', 'Interpreter', 'Latex');
        ylabel('$y / m$', 'Interpreter', 'Latex');
        title([field, ' / ', Unit, ' at ', ' $\omega_{LH}t = ', num2str(dt*(k-1)), '$'], 'Interpreter', 'Latex');
        set(gca, 'Fontsize', 12);
        savefig(Ffield, ['plot\', field, num2str(k), '.fig']);
        saveas(Ffield, ['plot\', field, num2str(k), '.png']);
        close(Ffield);
    end
end