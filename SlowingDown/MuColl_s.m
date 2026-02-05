function [mus, mupars, mupers, muepsis] = MuColl_s(Ph, Ps, v)
% Needs BasicPlasma/*.m
    v = v(:);
    e = 1.602176565e-19;
    epsilon_0 = 8.854187817e-12;
    if ~isfield(Ph, 'T')
        Ph.E = Eslowdown(Ph.m, Ph.u0, Ph.uc);
        Ph.T = 2*Ph.E/3/e;
    end
    mus = zeros(numel(v), numel(Ps(1).T), numel(Ps));
    mupars = zeros(numel(v), numel(Ps(1).T), numel(Ps));
    mupers = zeros(numel(v), numel(Ps(1).T), numel(Ps));
    muepsis = zeros(numel(v), numel(Ps(1).T), numel(Ps));
    for k = 1:numel(Ps)
        Ps(k).vt = sqrt(2*Ps(k).T*e/Ps(k).m);
        % PhPs = PhjPs(Ph, Ps(k));
        Ps(k).LnL = LnLambda(Ph, Ps(k));
        Ps(k).Gamma12 = Ps(k).n0.*(Ps(k).q*Ph.q*e^2/epsilon_0/Ph.m)^2/4/pi.*Ps(k).LnL;
        x1 = v./Ps(k).vt;
        x12 = x1.^2;
        [mux, muxp, ~] = muerr(x12);
        G2vt3 = Ps(k).vt.^(-3) .* Ps(k).Gamma12;
        mus(:, :, k) = G2vt3 .* (1 + Ph.m/Ps(k).m) .* mux.*x1.^(-3);
        mupars(:, :, k) = G2vt3 .* mux.*x1.^(-5);
        mupers(:, :, k) = 2*G2vt3 .* x1.^(-3) .* (mux + muxp - 0.5*mux.*x1.^(-2));
        muepsis(:, :, k) = 2*G2vt3 .* x1.^(-3) .* (Ph.m/Ps(k).m.*mux - muxp);
    end
end
% 
% function PhPs = PhjPs(Ph, Ps)
% % Convert a fast particle struct to a background particle struct.
%     S = fieldnames(Ps);
%     PhPs = Ps;
%     PhPs = [PhPs, PhPs(end)];
%     for k = 1:numel(S)
%         PhPs(end).(S{k}) = Ph.(S{k});
%     end
% end