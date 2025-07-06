function varargout = FsolveKperp(w, k_para, k_perp_ref, B, Ps, fun)
% Calculate k_perp with accurate Maxwellian dielectric tensor.
% k_perp_ref : reference k_perp
% fun : function handle, determine the dispersion relation to be used.
    f = @(k) fun(w, k_para, k, B, Ps);
    
    %----------------Calc with detailed steps----------------%
    % options = optimset('Display','iter');
    % [k_perp, ReD, exitflag, output] = fsolve(f, k_perp_ref, options);
    
    %----------------Calc without std output----------------%
    options = optimset('Display','none');
    [k_perp, D, exitflag, output] = fsolve(f, k_perp_ref*(1 + 0.0000001*1i), options);

    %----------------Calc with std output----------------%
    % [k_perp, D, exitflag, output] = fsolve(f, k_perp_ref*(1 + 0.0000001*1i));
    
    %----------------Send error message to diary----------------%
    if exitflag == 0
        fprintf('FsolveNperp : Iteration exceeds at k_perp_ref = %f\n', k_perp_ref);
    % elseif exitflag ~= -1
    %     fprintf('FzeroNperp : fzero failed at k_perp_ref = %f\n', k_perp_ref);
    end
    
    %----------------Output----------------%
    varargout{1} = k_perp;
    varargout{2} = exitflag;
    varargout{3} = D;
    varargout{4} = output;
end