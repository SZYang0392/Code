function varargout = setprof(Ps, r, Opts)
    % Create genray_profs_in.nc with profile parameters
        arguments
            % Ps is a struct array with fields : 
            % q(in e), m(in kg), n0(array, in m^-3), T(array, in eV)
            % Ps(1) should be electrons
            Ps
            % r : the transport radial mesh (?)
            r
            % qnorm : adjust density to keep charge neutrality
            Opts.qnorm = true
            % maxiter : maximum iteration number in charge neutralizing
            Opts.maxiter = 1000
            % Qerrmax : maximum relative error of charge
            Opts.Qerrmax = 1e-4
            % gfilename : the name of the g eqdsk file to be used in GENRAY
            Opts.gfilename = 'gfile_EFIT'
            % filename : the name of the generated .nc GENRAY input file.
            Opts.filename = 'genray_profs_in.nc'
            % Zeff : Zeff profile
            Opts.Zeff
        end
    
        %% Prepare the data
        % ---------------- Dimensions ---------------- %
        data.nspecgr = numel(Ps);
        data.nj = numel(Ps(1).n0);
        % ---------------- Create arrays ---------------- %
        data.charge = nan(data.nspecgr, 1);
        data.dmass = nan(data.nspecgr, 1);
        data.en = nan(data.nj, data.nspecgr);
        data.r = nan(data.nj, 1);
        data.temp = nan(data.nj, data.nspecgr);
        data.Zeff = nan(data.nj, 1);
        % ---------------- gfile name ---------------- %
        data.eqdsk_name = repmat(' ', 1, 256);
        data.eqdsk_name(1:numel(Opts.gfilename)) = Opts.gfilename;
        data.eqdsk_name = data.eqdsk_name.';
    
        % ---------------- Write profiles ---------------- %
        % r
        if isscalar(r)
            data.r = linspace(0, r, data.nj).';
        else
            data.r = r(:);
        end
        % Particle profiles
        for k = 1:data.nspecgr
            data.charge(k) = abs(Ps(k).q);
            data.dmass(k) = Ps(k).m./Ps(1).m;
            data.en(:, k) = Ps(k).n0(:).*1e-6;
            data.temp(:, k) = Ps(k).T(:).*1e-3;
        end
        if Opts.qnorm
            Qe = data.en(:, 1);
            Qi = data.en(:, 2:end)*data.charge(2:end);
            Qerr = abs(Qi./Qe - 1);
            niter = 0;
            while niter <= Opts.maxiter && ~all(Qerr <= Opts.Qerrmax, 'all')
                data.en(:, 2:end) = data.en(:, 2:end).*Qe./Qi;
                Qi = data.en(:, 2:end)*data.charge(2:end);
                Qerr = abs(Qi./Qe - 1);
                niter = niter + 1;
            end
            if niter > Opts.maxiter && ~all(Qerr <= Opts.Qerrmax, 'all')
                fprintf('===========================================================\n');
                fprintf('Error in GENRAY.setprof : Unable to keep charge neutrality.\n');
                fprintf('===========================================================\n');
            end
        end
        % Sort by charge and Temperature
        Ttemp = mean(data.temp(:, 2:end), 1).';
        [~, I] = sortrows([data.charge(2:end), Ttemp], [1, 2]);
        I = I + 1;
        data.charge(2:end) = data.charge(I);
        data.dmass(2:end) = data.dmass(I);
        data.en(:, 2:end) = data.en(:, I);
        data.temp(:, 2:end) = data.temp(:, I);
        % Zeff (2007_P.Phys.Fus.Ener_J.Freidberg, eq.3.41)
        if isfield(Opts, 'Zeff')
            data.Zeff = Opts.Zeff(:);
        else
            data.Zeff = data.en(:, 2:end)*(data.charge(2:end).^2)./data.en(:, 1);
        end
    
    
        %% Write the .nc file
        % ---------------- Create a new file ---------------- %
        ncid = netcdf.create(Opts.filename, 'CLOBBER');
    
        % ---------------- Create global attributes ---------------- %
        varid = netcdf.getConstant('GLOBAL');
        netcdf.putAtt(ncid, varid, 'title', 'Profile data passed from ONETWO to GENRAY');
    
        % ---------------- Create dimensions ---------------- %
        dimid = nan(1, 4);
        dimid(1) = netcdf.defDim(ncid, 'char256dim', 256);
        dimid(2) = netcdf.defDim(ncid, 'nj_dim', data.nj);
        dimid(3) = netcdf.defDim(ncid, 'nspecgr_dim', data.nspecgr);
        dimid(4) = netcdf.defDim(ncid, 'unity', 1);
    
        %% Create variables
        varid = nan(1, 9);
        % Define
        varid(1) = netcdf.defVar(ncid, 'charge', 'NC_DOUBLE', dimid(3));
        varid(2) = netcdf.defVar(ncid, 'dmass', 'NC_DOUBLE', dimid(3));
        varid(3) = netcdf.defVar(ncid, 'en', 'NC_DOUBLE', dimid([2, 3]));
        varid(4) = netcdf.defVar(ncid, 'eqdsk_name', 'NC_CHAR', dimid(1));
        varid(5) = netcdf.defVar(ncid, 'nj', 'NC_INT', dimid(4));
        varid(6) = netcdf.defVar(ncid, 'nspecgr', 'NC_INT', dimid(4));
        varid(7) = netcdf.defVar(ncid, 'r', 'NC_DOUBLE', dimid(2));
        varid(8) = netcdf.defVar(ncid, 'temp', 'NC_DOUBLE', dimid([2, 3]));
        varid(9) = netcdf.defVar(ncid, 'zeff', 'NC_DOUBLE', dimid(2));
        netcdf.endDef(ncid);
        % Write data
        netcdf.putVar(ncid, varid(1), data.charge);
        netcdf.putVar(ncid, varid(2), data.dmass);
        netcdf.putVar(ncid, varid(3), data.en);
        netcdf.putVar(ncid, varid(4), data.eqdsk_name);
        netcdf.putVar(ncid, varid(5), data.nj);
        netcdf.putVar(ncid, varid(6), data.nspecgr);
        netcdf.putVar(ncid, varid(7), data.r);
        netcdf.putVar(ncid, varid(8), data.temp);
        netcdf.putVar(ncid, varid(9), data.Zeff);
        % Put attributes
        netcdf.reDef(ncid);
        netcdf.putAtt(ncid, varid(1), 'long_name', 'Charge number of e (pos) and ion species');
        netcdf.putAtt(ncid, varid(2), 'long_name', 'Atomic wt of each spec, normd to elec mass');
        netcdf.putAtt(ncid, varid(3), 'long_name', 'Electron (first) and ion densities  ');
        netcdf.putAtt(ncid, varid(3), 'units', '/cm**3');
        netcdf.putAtt(ncid, varid(4), 'long_name', 'Name of input eqdsk; may vary with time');
        netcdf.putAtt(ncid, varid(5), 'long_name', 'Transport radial mesh dimension ');
        netcdf.putAtt(ncid, varid(6), 'long_name', 'Number of plasma species, incl electrons');
        netcdf.putAtt(ncid, varid(7), 'long_name', 'transport radial mesh');
        netcdf.putAtt(ncid, varid(7), 'units', 'cms');
        netcdf.putAtt(ncid, varid(8), 'long_name', 'Electron (first) and ion temperatures      ');
        netcdf.putAtt(ncid, varid(8), 'units', 'keV');
        netcdf.putAtt(ncid, varid(9), 'long_name', 'zeff on transport radial mesh');
        netcdf.endDef(ncid);
    
        %% Close file
        netcdf.close(ncid);
    end