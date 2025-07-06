function varargout = Polar(w, k_para, k_perp, B, Ps, type, tol)
    % Calculate polarization
    % B, w, N_para, N_perp are scalars
    % f : function handle of the susceptibility tensor
    % E : [perp; 0; para]

    % For Maxwellian initial distribution only.
    f = @KhiMaxwell;

    c = 299792458;
    N_para = k_para*c/w;
    N_perp = k_perp*c/w;
    N = [N_perp, 0, N_para];

    % Calculate susceptibility tensor
    Khis = zeros(3, 3, numel(Ps));
    for k = 1:numel(Ps)
        Khis(:, :, k) = f(w, k_para, k_perp, B, Ps(k), 0);
    end
    Khi = sum(Khis, 3);

    % Calculate dielectric tensor
    if strcmpi(type, 'es')         % Compare strings, capitalization ignored.
        D = Khi + eye(3);
    elseif strcmpi(type, 'em')
        D = N.'*N - eye(3).*(N_para.^2 + N_perp.^2) + eye(3) + Khi;
    else                           % Both ES and EM
        D(:, :, 1) = Khi + eye(3);
        D(:, :, 2) = N.'*N - eye(3).*(N_para.^2 + N_perp.^2) + Khi + eye(3);
    end

    % Calculate polarization
    if all(~isnan(D), 'all') && all(~isinf(D), 'all')
        if strcmpi(type, 'es') || strcmpi(type, 'em')
            if nargin >= 7
                E = null(D, tol);
            else
                E = null(D);
            end
        else                           % Both ES and EM
            if nargin >= 7
                E(:, :, 1) = null(D(:, :, 1), tol);
                E(:, :, 2) = null(D(:, :, 2), tol);
            else
                E(:, :, 1) = null(D(:, :, 1));
                E(:, :, 2) = null(D(:, :, 2));
            end
        end
    else
        E = [];
    end


    varargout{1} = E;
    varargout{2} = D;
    varargout{3} = Khi;
    varargout{4} = Khis;
end