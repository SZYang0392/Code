function data = load_flux(filename)
%CQL.load_flux  Load CQL3D distribution function from a _flux_*.nc file.
%
%   data = CQL.load_flux(filename) reads a CQL3D _flux_N.nc file (written
%   when netcdfshort=''enabled'').  The distribution function f is
%   converted to SI units.
%
%   The flux files contain the full distribution function
%   f(ydimf, xdimf, rdim) at a single saved time step.
%
%   Example:
%       f1 = CQL.load_flux('CFETR_NB_transport_flux_1.nc');
%       imagesc(f1.x, f1.y, log10(abs(squeeze(f1.f(:, :, end))')));

    ncid = netcdf.open(filename, 'NC_NOWRITE');
    [~, numvars, ~, ~] = netcdf.inq(ncid);

    data = struct();
    data.filename = filename;                     %#ok<STRNU>

    for i = 1:numvars
        varid = i - 1;
        val   = netcdf.getVar(ncid, varid);
        [vname, ~, ~, ~] = netcdf.inqVar(ncid, varid);
        data.(vname) = val;
    end

    netcdf.close(ncid);

    % Convert vnorm (cm/s -> m/s)
    if isfield(data, 'vnorm')
        data.vnorm = data.vnorm .* 0.01;          % cm/s -> m/s
    end

    % Convert rhomax (cm -> m)
    data.rhomax = data.rhomax .* 0.01;          % cm -> m

    % Convert f and favr_thet0 from cm^-3 to m^-3.
    % (The stored f is a number density: f_code = f_cgs * vnorm^3.)
    scale_f = 1e6;
    if isfield(data, 'f')
        data.f = data.f .* scale_f;
    end
    if isfield(data, 'favr_thet0')
        data.favr_thet0 = data.favr_thet0 .* scale_f;
    end

    clear ncid numvars i varid val vname scale_f;
end