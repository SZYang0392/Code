function Fk = fftkfield(F, Opts)
% Calculate the Forier transform of sepcific fields.
    arguments
        F 
        Opts.spec string = 'all'
        Opts.nx (1, 2) {mustBeInteger, mustBeNonnegative}
        Opts.ny (1, 2) {mustBeInteger, mustBeNonnegative}
        Opts.nz (1, 2) {mustBeInteger, mustBeNonnegative}
    end

    % Select domain
    Fname = fieldnames(F);
    if isfield(Opts, 'nx')
        F.nx = Opts.nx(end) - Opts.nx(1) + 1;
        F.x = F.x(Opts.nx(1):Opts.nx(end));
        for j = 14:numel(Fname)
            F.(Fname{j}) = F.(Fname{j})(Opts.nx(1):Opts.nx(end), :, :);
        end
    end
    if isfield(Opts, 'ny')
        F.ny = Opts.ny(end) - Opts.ny(1) + 1;
        F.y = F.y(Opts.ny(1):Opts.ny(end));
        for j = 14:numel(Fname)
            F.(Fname{j}) = F.(Fname{j})(:, Opts.ny(1):Opts.ny(end), :);
        end
    end
    if isfield(Opts, 'nz')
        F.nz = Opts.nz(end) - Opts.nz(1) + 1;
        F.z = F.z(Opts.nz(1):Opts.nz(end));
        for j = 14:numel(Fname)
            F.(Fname{j}) = F.(Fname{j})(:, :, Opts.nz(1):Opts.nz(end));
        end
    end
    % Parameters of each dim.
    fsize = [F.nx, F.ny, F.nz];
    dimexist = fsize > 1;
    dx = mean(F.x(2:end) - F.x(1:end-1), 'all');
    dy = mean(F.y(2:end) - F.y(1:end-1), 'all');
    dz = mean(F.z(2:end) - F.z(1:end-1), 'all');
    dr = [dx, dy, dz];

    % Calculate frequency domain
    r = {F.x, F.y, F.z};
    DK = zeros(1, 3);
    Kdomain = cell(1, 3);
    for k = 1:3
        if dimexist(k)
            dkk = 2*pi/fsize(k)/dr(k);
            DK(k) = dkk;
            Kdomain{k} = ((1:fsize(k)) - (fsize(k) + 1)/2.0)*dkk;
            Kdomain{k} = reshape(Kdomain{k}, size(r{k}));
        else
            Kdomain{k} = 0;
        end
    end

    % Perform Fourier transform
    iffft = 0;
    fn = fieldnames(F);
    nfn = numel(fn);
    if strcmp(Opts.spec, 'all')
        for k = nfn:-1:nfn-F.nfield+1
            FieldX = F.(fn{k});
            FieldK = fftcalc(FieldX, dr);
            F = rmfield(F, fn{k});
            F.([fn{k}, 'k']) = FieldK;
            iffft = 1;
        end
    else
        for k = nfn:-1:nfn-F.nfield+1
            if strcmp(Opts.spec, fn{k})
                FieldX = F.(fn{k});
                FieldK = fftcalc(FieldX, dr);
                F = rmfield(F, fn{k});
                F.([fn{k}, 'k']) = FieldK;
                iffft = 1;
            end
        end
    end

    % Check if Fourier transform is done.
    if ~iffft
        fprintf('fft not executed \n');
    end

    % Set Fourier field Fk.
    Fk = F;
    Fk.dr = dr;
    Fk.dx = dx;
    Fk.dy = dy;
    Fk.dz = dz;
    Fk.dk = DK;
    Fk.dkx = DK(1);
    Fk.dky = DK(2);
    Fk.dkz = DK(3);
    Fk.kdomain = Kdomain;
    Fk.kx = Kdomain{1};
    Fk.ky = Kdomain{2};
    Fk.kz = Kdomain{3};
end

function k = fftcalc(x, dr)
% Multidimension Fourier tramsform.
% x : Field to be transformed.
% dr : Intervals of each dimension.
    k = x;
    for j = 1:3
        sizex = size(k, j);
        if sizex > 1
            k = fft(k, sizex, j)*dr(j);
        end
    end
    k = fftshift(k);
end