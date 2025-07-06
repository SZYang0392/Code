function varargout = KESNM(w, k_para, k_perp, P)
    %--------------Physical constants--------------%
    epsilon_0 = 8.854187817e-12;
    e = 1.602176565e-19;
    %--------------Check Drift--------------%
    if isfield(P, 'v0')
        P.v0 = 0;
    end
    %--------------Susceptibility--------------%
    P.wp2 = (P.q.*e).^2.*P.n0./P.m./epsilon_0;
    if ~isfield(P, 'v0')
        P.v0 = 0;
    end
    w1 = w - k_para.*P.v0;
    if isfield(P, 'uniso')
        if P.uniso == 1
            if isfield(P, 'T_para')
                P.v_para = sqrt(2*e.*P.T_para./P.m);
            else
                P.v_para = sqrt(2*e.*P.T./P.m);
            end
            if isfield(P, 'T_perp')
                P.v_perp = sqrt(2*e.*P.T_perp./P.m);
            else
                P.v_perp = sqrt(2*e.*P.T./P.m);
            end
            kvt = sqrt(k_para.^2.*P.v_para.^2 + k_perp.^2.*P.v_perp.^2);
        else
            P.vt = sqrt(2*e.*P.T./P.m);
            P.kvt = P.vt.*sqrt(k_para.^2 + k_perp.^2);
        end
    else
        P.vt = sqrt(2*e.*P.T./P.m);
        P.kvt = P.vt.*sqrt(k_para.^2 + k_perp.^2);
    end
    P.zeta = w1./P.kvt;
    P.Khi = 2*P.wp2./(w1.^2).*Gesnm(P.zeta);
    varargout{1} = P.Khi;
    varargout{2} = P;
end

function G = Gesnm(z)
    G = zeros(size(z));
    ind = isinf(z);
    G(ind) = -0.5;
    G(~ind) = z(~ind).^2.*(1 + z(~ind).*Z_plasma(z(~ind)));
end