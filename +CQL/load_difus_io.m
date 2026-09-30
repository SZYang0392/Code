function data = load_difus_io(filename)
%CQL.load_difus_io  Load CQL3D radial diffusion I/O netCDF (_difus_io.nc).
%
%   data = CQL.load_difus_io(filename) reads a CQL3D _difus_io.nc file
%   and returns a struct with variables converted to SI units.
%
%   filename can be either:
%     - The _difus_io.nc file directly, or
%     - The main mnemonic.nc file (the _difus_io suffix is appended).
%
%   Key variables:
%     d_rr  -- radial diffusion coefficient            cm^2/sec -> m^2/s
%     d_r   -- radial pinch velocity (if present)       cm/sec   -> m/s
%     rya   -- normalised radial mesh at bin centres
%     x     -- normalised momentum-per-mass grid
%     y     -- pitch angle grid (radians)
%
%   Example:
%       d = CQL.load_difus_io('CFETR_NB_transport.nc');
%       mesh(d.x, d.rya, squeeze(d.d_rr(1,:,:))');  % d_rr at y=1 vs x, rya

    % Derive _difus_io.nc name if a main mnemonic.nc is given
    [p, name, ext] = fileparts(filename);
    if ~contains(name, '_difus_io')
        filename = fullfile(p, [name '_difus_io' ext]);
    end

    % -----------------------------------------------------------------
    % 1.  Open netCDF, read every variable, and collect unit metadata
    % -----------------------------------------------------------------
    ncid = netcdf.open(filename, 'NC_NOWRITE');
    [~, numvars, ~, ~] = netcdf.inq(ncid);

    data = struct();
    data.filename = filename;                     %#ok<STRNU>

    var_units = cell(0, 2);

    for i = 1:numvars
        varid = i - 1;
        val   = netcdf.getVar(ncid, varid);
        [vname, ~, ~, ~] = netcdf.inqVar(ncid, varid);
        data.(vname) = val;
        try
            u = netcdf.getAtt(ncid, varid, 'units');
            if ischar(u) || isstring(u)
                var_units(end+1, :) = {vname, char(u)}; %#ok<AGROW>
            end
        catch
        end
    end

    netcdf.close(ncid);

    % -----------------------------------------------------------------
    % 2.  Convert CGS -> SI
    % -----------------------------------------------------------------
    % Reuse the same conversion engine from the main load, but add
    % difus_io-specific unit strings.
    data = apply_unit_conversions_difus(data, var_units);
    % Convert rpconz
    data.rpconz = data.rpconz .* 0.01;

    clear ncid numvars i varid val vname u;
end

% =========================================================================
function data = apply_unit_conversions_difus(data, var_units)
% Same logic as apply_unit_conversions, with difus_io-specific additions.

    persistent lut;
    if isempty(lut)
        lut = { ...
            % ---- Length ----
            'cms',                      0.01; ...
            'cms^2',                    1e-4; ...
            'cms^3',                    1e-6; ...
            'cms/sec',                  0.01; ...
            ... % ---- Diffusion / pinch ----
            'cm**2/se',                 1e-4; ...  % cm^2/s -> m^2/s (truncated)
            'cm**2/sec',                1e-4; ...  % cm^2/s -> m^2/s
            'cm/sec',                   0.01; ...  % cm/s   -> m/s
            ... % ---- Energy ----
            'kev',                      1e3; ...   % keV    -> eV
            };
    end

    si_or_dimless = {'', 'none', 'unitless', 'radians', 'seconds', 'secs', ...
                     'amps', 'watts', 'norm'};

    for k = 1:size(var_units, 1)
        vname    = var_units{k, 1};
        unit_str = strtrim(lower(var_units{k, 2}));

        if ~isfield(data, vname)
            continue;
        end
        val = data.(vname);
        if ~isnumeric(val) && ~islogical(val)
            continue;
        end
        if any(strcmp(unit_str, si_or_dimless))
            continue;
        end

        factor = 1;
        for r = 1:size(lut, 1)
            if strcmp(unit_str, lut{r, 1})
                factor = lut{r, 2};
                break;
            end
        end
        if factor == 1
            for r = 1:size(lut, 1)
                pattern = lut{r, 1};
                if length(unit_str) >= length(pattern) ...
                   && strncmp(unit_str, pattern, length(pattern))
                    factor = lut{r, 2};
                    break;
                end
            end
        end

        if factor ~= 1
            data.(vname) = double(val) .* factor;
        end
    end
end
