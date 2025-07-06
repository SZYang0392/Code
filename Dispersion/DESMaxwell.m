function varargout = DESMaxwell(w, k_para, k_perp, B, Ps, ifstrict)
% Calculate ES Dispersion relation with Maxwellian distribution.
% w, N_para, N_perp, B, Ps(k).* are scalars.

    % Default input parameters
    if nargin < 6
        ifstrict = false;
    end

    c = 299792458;
    N_para = k_para.*c./w;
    N_perp = k_perp.*c./w;
    N = [N_perp, 0, N_para];
    Khis = zeros(3, 3, numel(Ps));
    Ds = zeros(numel(Ps), 1);
    if ifstrict
        %% ====================== (Chi + I)*N = 0 ====================== %%
        for k = 1:numel(Ps)
            if isfield(Ps(k), 'mag')
                if Ps(k).mag == 0
                    Khis(:, :, k) = KESNM(w, k_para, k_perp, Ps(k))*eye(3);
                else
                    Khis(:, :, k) = KhiMaxwell(w, k_para, k_perp, B, Ps(k), 0);
                end
            else
                Khis(:, :, k) = KhiMaxwell(w, k_para, k_perp, B, Ps(k), 0);
            end
            Ds(k) = N*Khis(:, :, k)*(N.');
        end
        Khi = sum(Khis, 3);
        Dielectric = Khi + eye(3);
        D0 = Dielectric*(N.');
        varargout{1} = (D0')*D0;
    else
        %% ====================== N*(Chi + I)*N = 0 ====================== %%
        for k = 1:numel(Ps)
            if isfield(Ps(k), 'mag')
                if Ps(k).mag == 0
                    Khis(:, :, k) = KESNM(w, k_para, k_perp, Ps(k))*eye(3);
                else
                    Khis(:, :, k) = KhiMaxwell(w, k_para, k_perp, B, Ps(k), 0);
                end
            else
                Khis(:, :, k) = KhiMaxwell(w, k_para, k_perp, B, Ps(k), 0);
            end
            Ds(k) = N*Khis(:, :, k)*(N.');
        end
        varargout{1} = sum(Ds, 'all') + N*eye(3)*(N.');
        Khi = sum(Khis, 3);
    end
    varargout{2} = Ds;
    varargout{3} = Khi;
    varargout{4} = Khis;
    %--------------Check Units--------------%
    N2 = N_para.^2 + N_perp.^2;
    varargout{1} = varargout{1}.*N2;
    varargout{2} = varargout{2}.*N2;
end