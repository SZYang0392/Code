function Khi = KhiESMaxwell(w, k_para, k_perp, P)
    % Unmagnetizd electromagnetic dielectric tensor
    %--------------Physical constants--------------%
    epsilon_0 = 8.854187817e-12;
    e = 1.602176565e-19;
    %--------------Calculate parameters--------------%
    P.wp2 = (P.q.*e).^2.*P.n0./P.m./epsilon_0;
    if ~isfield(P, 'v0')
        P.v0 = 0;
    end
    w1 = w - k_para.*P.v0;
    % kvt
    if ~isfield(P, 'uniso')
        P.T_para = P.T;
        P.T_perp = P.T;
    end
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
    kvt = sqrt(k_para.^2.*P.v_para.^2 + k_perp.*P.v_perp.^2);
    ktot = sqrt(k_para.^2 + k_perp.^2);
    vt = kvt./ktot;
    % zeta
    z = w1./kvt;
    %--------------Martices--------------%
    % Use Normalized parameters
    A = diag([P.v_perp, P.v_perp, P.v_para])./vt;
    P.v0 = P.v0./vt;

    coskv = k_perp.*P.v_perp./kvt;
    sinkv = k_para.*P.v_para./kvt;
    B = [coskv, 0, sinkv ;...
         0,     1, 0     ;...
         sinkv, 0, -coskv];
    AB = A*B;
    K1 = AB*diag([Z3(z), 0.5*Z1(z), 0.5*Z1(z)])*(AB.');
    K2 = AB*[0, 0, Z2(z)*P.v0;...
             0, 0, 0         ;...
             0, 0, 0         ];
    K2 = K2 + K2.';
    K3 = [0, 0, 0;...
          0, 0, 0;...
          0, 0, P.v0.^2.*Z1(z)];
    param = 0.5*k_para.*k_perp./(ktot.^2).*(P.v_para./P.v_perp - P.v_perp./P.v_para);
    K4 = AB*[0,            0, Z1(z).*param;...
             0,            0, 0;...
             Z1(z).*param, 0, 0]*(AB.');
    K5 = AB*[0, 0, 0;...
             0, 0, 0;...
             0, 0, P.v0.*param];
    K5 = K5 + K5.';
    Khi = P.wp2./(w.^2).*(-eye(3) + 2.*(K1 + K2 + K3 + K4 + K5));
end

function Z = Z3(z)
    Z = 0.5 + z.^2.*Z1(z);
end
function Z = Z2(z)
    Z = z.*Z1(z);
end
function Z = Z1(z)
    Z = 1 + z.*Z_plasma(z);
end