function plotE(E, field, tscale, mime)
% field : cell array with name of struct fields.
% mime : mi/me.
    if nargin < 4
        mime = 1836;
    end
    if nargin < 3
        tscale = 1;
        %tscale = w/wce
        t = E.stime;
    else
        t = E.stime*tscale*mime;
    end

    Ffield = figure;
    hold on;
    for k = 1:numel(field)
        F = E.(field{k});
        plot(t, F/F(1), 'LineWidth', 1.5, 'DisplayName', field{k});
    end
    legend;
    if nargin >= 3
        xlabel('$\omega_{LH}t$', 'Interpreter', 'Latex');
    else
        xlabel('$\omega_{ci}t$', 'Interpreter', 'Latex');
    end
    ylabel('$E/E(0)$', 'interpreter', 'Latex');
    savefig(Ffield, [field{1}, '.fig']);
    saveas(Ffield, [field{1}, '.png']);
    % close(Ffield);
end