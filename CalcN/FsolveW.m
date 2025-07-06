function varargout = FsolveW(w_ref, k_para, k_perp, B, Ps, fun)
% Calculate k_perp with accurate Maxwellian dielectric tensor.
% w : reference w
% fun : function handle, determine the dispersion relation to be used.
    f = @(w) fun(w, k_para, k_perp, B, Ps);
    
    %----------------Calc with detailed steps----------------%
    % options = optimset('Display','iter');
    % [w, D, exitflag, output] = fsolve(f, w_ref, options);
    
    %----------------Calc without std output----------------%
    options = optimset('Display','none');
    [w, D, exitflag, output] = fsolve(f, w_ref, options);

    %----------------Calc with std output----------------%
    % [w, D, exitflag, output] = fsolve(f, w_ref);
    
    %----------------Send error message to diary----------------%
    if exitflag == 0
        fprintf('FsolveNperp : Iteration exceeds at w_ref = %f\n', w_ref);
    % elseif exitflag ~= -1
    %     fprintf('FzeroNperp : fzero failed at k_perp = %f\n', k_perp);
    end
    
    %----------------Output----------------%
    varargout{1} = w;
    varargout{2} = exitflag;
    varargout{3} = D;
    varargout{4} = output;
end