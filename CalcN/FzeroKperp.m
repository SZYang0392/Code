function varargout = FzeroKperp(w, k_para, k_perp_ref, B, Ps, fun)
% Calculate k_perp with accurate Maxwellian dielectric tensor.
% k_perp_ref : reference k_perp
% fun : function handle, determine the dispersion relation to be used.
    f = @(k) real(fun(w, k_para, k, B, Ps));

    %----------------Calc with detailed steps----------------%
    % options = optimset('Display','iter');
    % [k_perp, ReD, exitflag, output] = fzero(f, k_perp_ref, options);

    %----------------Calc without detailed steps----------------%
    options = optimset('Display','none');
    [k_perp, ReD, exitflag, output] = fzero(f, k_perp_ref, options);

    %----------------Send error message to diary----------------%
    if exitflag == -3
        fprintf('FzeroNperp : NaN or Inf occur at k_perp_ref = %f\n', k_perp_ref);
    elseif exitflag == -4
        fprintf('FzeroNperp : Complex values occur at k_perp_ref = %f\n', k_perp_ref);
    elseif exitflag == -6
        fprintf('FzeroNperp : No sign change detected around k_perp_ref = %f\n', k_perp_ref);
    % elseif exitflag ~= -1
    %     fprintf('FzeroNperp : fzero failed at k_perp_ref = %f\n', k_perp_ref);
    end

    %----------------Output----------------%
    varargout{1} = k_perp;
    varargout{2} = exitflag;
    varargout{3} = ReD;
end