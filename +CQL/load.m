function data = load(filename)
%CQL.load  Read CQL3D netCDF output file and convert to SI units.
%
%   data = CQL.load(filename) reads a CQL3D main netCDF output file and
%   returns a struct containing all variables.  Units are converted from
%   the native CQL3D mix of CGS/SI to SI: lengths in m, magnetic field in
%   T, energy/temperature in eV, current density in A/m^2, power density
%   in W/m^3, poloidal flux in Wb, resistivity in ohm*m, etc.
%
%   Example:
%       d = CQL.load('CFETR_NB_LHW.nc');
%       plot(d.rya, d.temp(:, end));  % Temperature profile at last time step

    % ---------------------------------------------------------------------
    % 1.  Open netCDF, read every variable, and collect unit metadata
    % ---------------------------------------------------------------------
    ncid = netcdf.open(filename, 'NC_NOWRITE');
    [~, numvars, ~, ~] = netcdf.inq(ncid);

    data = struct();
    data.filename = filename;                     %#ok<STRNU>

    % var_units{n,1} = variable name, var_units{n,2} = units string
    var_units = cell(0, 2);

    for i = 1:numvars
        varid = i - 1;                           % netCDF uses 0-based IDs
        val   = netcdf.getVar(ncid, varid);
        [vname, ~, ~, ~] = netcdf.inqVar(ncid, varid);
        data.(vname) = val;

        % Attempt to read the "units" attribute
        try
            u = netcdf.getAtt(ncid, varid, 'units');
            if ischar(u) || isstring(u)
                var_units(end+1, :) = {vname, char(u)}; %#ok<AGROW>
            end
        catch
        end
    end

    netcdf.close(ncid);

    % ---------------------------------------------------------------------
    % 2.  Apply CGS / mixed-unit  ->  SI conversions
    % ---------------------------------------------------------------------
    data = apply_unit_conversions(data, var_units);

    % ---------------------------------------------------------------------
    % 3.  Distribution-function rescaling (cm^-3 -> m^-3)
    % ---------------------------------------------------------------------
    % The stored f and favr_thet0 are number densities in cm^-3
    % (diaggnde.f: f_code/vnorm^3 = f_cgs, so f_code has units cm^-3).
    % Converting to SI therefore only requires cm^-3 -> m^-3.
    scale_f = 1e6;
    if isfield(data, 'f')
        data.f = data.f .* scale_f;
    end
    if isfield(data, 'favr_thet0')
        data.favr_thet0 = data.favr_thet0 .* scale_f;
    end

    % ---------------------------------------------------------------------
    % 4.  Fix compound-unit variables (e.g. rfpwr has both W/cm^3 and W)
    % ---------------------------------------------------------------------
    data = fix_compound_units(data);

    clear ncid numvars i varid val vname u scale_f;
end

% =========================================================================
function data = apply_unit_conversions(data, var_units)
% Apply multiplicative conversion factors based on the "units" attribute.
% Handles exact matches and prefix matches for compound unit strings.

    % Look-up table: unit-pattern  ->  conversion factor  (-> SI)
    persistent lut;
    if isempty(lut)
        lut = { ...
            % ---- Length ----
            'cms',                      0.01; ...   % cm      -> m
            'cms^2',                    1e-4; ...   % cm^2    -> m^2
            'cms^3',                    1e-6; ...   % cm^3    -> m^3
            'cms/sec',                  0.01; ...   % cm/s    -> m/s
            ... % ---- Magnetic field ----
            'gauss',                    1e-4; ...   % G       -> T
            ... % ---- Flux (poloidal, CGS => Mx = G*cm^2 => 1 Mx = 1e-8 Wb) ----
            'cgs',                      1e-8; ...   % Mx      -> Wb
            ... % ---- Mass ----
            'grams',                    1e-3; ...   % g       -> kg
            'gram',                     1e-3; ...   % g       -> kg
            ... % ---- Current density ----
            'amps/cm**2',               1e4; ...    % A/cm^2  -> A/m^2
            'amps/cm^2',                1e4; ...
            'a/cm**2',                  1e4; ...
            'amps/cm**',                1e4; ...    % truncated 'amps/cm**2'
            ... % ---- Power density ----
            'watts/cm**3',              1e6; ...    % W/cm^3  -> W/m^3
            'watts/cm^3',               1e6; ...
            'w/cm^3',                   1e6; ...
            ... % ---- Number density ----
            '/cm**3',                   1e6; ...    % /cm^3   -> /m^3
            ... % ---- Electric field ----
            'volts/cm',                 1e2; ...    % V/cm    -> V/m
            ... % ---- Energy / temperature ----
            'kev',                      1e3; ...    % keV     -> eV
            ... % ---- Resistivity (CGS-Gaussian: resistivity has dimension of time) ----
            'cgs, seconds',             8.987551787e9; ... % s(CGS) -> ohm*m
            ... % ---- Compound strings (currv, pwrrf include integral notes) ----
            'amps/cm^2 (int:0,1 over dx =current density)', 1e4; ...
            'w/cm^3 (int:0,1 over dx =rf power density)',   1e6; ...
            'watts/cm**3, except watts for sorpwti',        1e6; ...  % rfpwr (bulk)
            };
    end

    % Unit strings that are already SI or dimensionless (skip conversion)
    si_or_dimless = {'', 'none', 'unitless', 'radians', 'seconds', 'secs', ...
                     'amps', 'watts', 'norm', 'normalized to vnorm', ...
                     'vnorm**3/(cm**3*(cm/sec)**3)'};

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

        % Skip if already SI / dimensionless
        if any(strcmp(unit_str, si_or_dimless))
            continue;
        end

        % Try exact-match look-up table
        factor = 1;
        for r = 1:size(lut, 1)
            if strcmp(unit_str, lut{r, 1})
                factor = lut{r, 2};
                break;
            end
        end

        % If no exact match, try prefix match (for compound strings)
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

% =========================================================================
function data = fix_compound_units(data)
% Fix variables whose "units" attribute describes a mix of units.
%
%   rfpwr:  columns 1:mrfn+2 are W/cm^3, column mrfn+3 is Watts.
%           After the bulk conversion they were all multiplied by 1e6,
%           so the last column must be divided back.

    % --- rfpwr -----------------------------------------------------------
    if isfield(data, 'rfpwr') && isfield(data, 'mrfn')
        mrfn = data.mrfn;
        ncol = size(data.rfpwr, 2);
        if mrfn > 0 && ncol >= mrfn + 3
            % Column mrfn+3 is integrated power (Watts), undo 1e6 factor
            data.rfpwr(:, mrfn+3) = data.rfpwr(:, mrfn+3) ./ 1e6;
        end
    end
end