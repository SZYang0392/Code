function data = load(filename, i_salphal)
    ncid = netcdf.open(filename);
    [~, numvars, ~, ~] = netcdf.inq(ncid);
    data.name = filename;
    for i=1:numvars
        temp = netcdf.getVar(ncid,i-1);
        [name,~,~,~]=netcdf.inqVar(ncid,i-1);
        data.(name) = temp;
    end
    netcdf.close(ncid);
    clear("ncid", "numvars", "temp", "name", "i");

    % Clear ineffective data
    Nray = numel(data.nrayelt);
    for k = 1:Nray
        Nrayelt = data.nrayelt(k);
        if Nrayelt >= size(data.ws, 1)
            continue;
        end
        data.ws(Nrayelt:end, k) = nan;
        data.seikon(Nrayelt:end, k) = nan;
        data.spsi(Nrayelt:end, k) = nan;
        data.wr(Nrayelt:end, k) = nan;
        data.wphi(Nrayelt:end, k) = nan;
        data.wz(Nrayelt:end, k) = nan;
        data.w_theta_pol(Nrayelt:end, k) = nan;
        data.wnpar(Nrayelt:end, k) = nan;
        data.wnper(Nrayelt:end, k) = nan;
        data.delpwr(Nrayelt:end, k) = 0;
        data.sdpwr(Nrayelt:end, k) =0;
        data.wdnpar(Nrayelt:end, k) = nan;
        data.cwexde(Nrayelt:end, k, :) = nan;
        data.cweyde(Nrayelt:end, k, :) = nan;
        data.cwezde(Nrayelt:end, k, :) = nan;
        data.fluxn(Nrayelt:end, k) = nan;
        data.sbtot(Nrayelt:end, k) = nan;
        data.sene(Nrayelt:end, k) = nan;
        data.ste(Nrayelt:end, k) = nan;
        data.szeff(Nrayelt:end, k) = nan;
        data.salphac(Nrayelt:end, k) = 0;
        data.salphac = abs(data.salphac);
        data.salphal(Nrayelt:end, k) = 0;
        data.salphal = abs(data.salphal);
        data.sb_r(Nrayelt:end, k) = nan;
        data.sb_z(Nrayelt:end, k) = nan;
        data.sb_phi(Nrayelt:end, k) = nan;
        data.wn_r(Nrayelt:end, k) = nan;
        data.wn_z(Nrayelt:end, k) = nan;
        data.wn_phi(Nrayelt:end, k) = nan;
        data.vgr_r(Nrayelt:end, k) = nan;
        data.vgr_z(Nrayelt:end, k) = nan;
        data.vgr_phi(Nrayelt:end, k) = nan;
        data.flux_z(Nrayelt:end, k) = nan;
        data.flux_r(Nrayelt:end, k) = nan;
        data.flux_phi(Nrayelt:end, k) = nan;
        data.w_eff_nc(Nrayelt:end, k) = nan;
        data.cweps11(Nrayelt:end, k, :) = nan;
        data.cweps12(Nrayelt:end, k, :) = nan;
        data.cweps13(Nrayelt:end, k, :) = nan;
        data.cweps21(Nrayelt:end, k, :) = nan;
        data.cweps22(Nrayelt:end, k, :) = nan;
        data.cweps23(Nrayelt:end, k, :) = nan;
        data.cweps31(Nrayelt:end, k, :) = nan;
        data.cweps32(Nrayelt:end, k, :) = nan;
        data.cweps33(Nrayelt:end, k, :) = nan;
    end
    if isfield(data, 'salphas')
        for k = 1:Nray
            Nrayelt = data.nrayelt(k) + 1;
            data.salphas(Nrayelt:end, k, :) = 0;
        end
        data.salphas = abs(data.salphas);
    end

    % Change to SI units
    c = 299792458;
    data.ws = data.ws/100;
    data.wr = data.wr/100;
    data.wz = data.wz/100;
    data.w_theta_pol = data.w_theta_pol*pi/180;
    data.delpwr = data.delpwr/1e7;
    data.deposite = [(data.delpwr(1, :)-data.delpwr(2, :))/2; ...
        (data.delpwr(1:end-2, :)-data.delpwr(3:end, :))/2; ...
        (data.delpwr(end-1, :)-data.delpwr(end, :))/2];
    % data.deposite = [zeros(1, data.nray); ...
    %     (data.delpwr(1:end-2, :)-data.delpwr(3:end, :))/2; zeros(1, data.nray)];
    data.sdpwr = data.sdpwr*100;
    data.fluxn = data.fluxn/(2.99792458^2*1e11);
    data.sbtot = data.sbtot/1e4;
    data.sene = data.sene*1e6;
    data.ste = data.ste*1e3; % in eV
    data.salphac = data.salphac*100;
    data.salphal = data.salphal*100;
    data.sb_r = data.sb_r/1e4;
    data.sb_z = data.sb_z/1e4;
    data.sb_phi = data.sb_phi/1e4;
    data.vgr_r = data.vgr_r*c;
    data.vgr_z = data.vgr_z*c;
    data.vgr_phi = data.vgr_phi*c;
    data.flux_r = data.flux_r*c;
    data.flux_z = data.flux_z*c;
    data.flux_phi = data.flux_phi*c;
    data.w_eff_nc = data.w_eff_nc*1e5;
    % data.dmas = data.dmas; % Normalized to electron mass
    data.w_tot_pow_absorb_at_refl_nc = data.w_tot_pow_absorb_at_refl_nc/1e7;
    data.power_inj_total = data.power_inj_total/1e7;
    data.power_total = data.power_total/1e7;
    data.powtot_e = data.powtot_e/1e7;
    data.powtot_i = data.powtot_i/1e7;
    data.powtot_cl = data.powtot_cl/1e7;
    data.voltot = data.voltot/1e6;
    data.areatot = data.areatot/1e4;
    data.pollentot = data.pollentot/100;
    data.binvol = data.binvol/1e6;
    data.binarea = data.binarea/1e4;
    data.pollen = data.pollen/100;
    data.densprof = data.densprof*1e6;
    data.temprof = data.temprof*1e3; % in eV
    data.spower = data.spower/1e7;
    data.w_dens_vs_r_nc = data.w_dens_vs_r_nc*1e6;
    data.w_temp_vs_r_nc = data.w_temp_vs_r_nc*1e3; % in eV
    data.s_cur_den_parallel = data.s_cur_den_parallel*1e4;
    data.s_cur_den_onetwo = data.s_cur_den_onetwo*1e4;
    data.s_cur_den_toroidal = data.s_cur_den_toroidal*1e4;
    data.s_cur_den_poloidal = data.s_cur_den_poloidal*1e4;
    data.powden = data.powden/10;
    data.powden_e = data.powden_e/10;
    data.powden_cl = data.powden_cl/10;
    data.powden_i = data.powden_i/10;
    if isfield(data, 'salphas')
        data.salphas = data.salphas*100;
    end
    if isfield(data, 'powtot_s')
        data.powtot_s = data.powtot_s/1e7;
    end
    if isfield(data, 'powden_s')
        data.powden_s = data.powden_s/10;
        temp = zeros(numel(data.rho_bin), size(data.powden_s, 2));
        temp(1:numel(data.powden_s)) = data.powden_s(:);
        data.powden_s = temp(1:end-1, :);
        data.powden_s = [data.powden_e, data.powden_s];
    end

    % Caculate nrho, mpol, ntol, krho, kpol, ktol
    data.wnrho = nan(size(data.spsi));
    data.mpol = nan(size(data.spsi));
    data.ntol = nan(size(data.spsi));
    data.krho = nan(size(data.spsi));
    data.kpol = nan(size(data.spsi));
    data.ktol = nan(size(data.spsi));
    w = 2*pi*data.freqcy;
    kr = data.wn_r*(w/c);
    kz = data.wn_z*(w/c);
    kphi = data.wn_phi*(w/c);
    bpol = sqrt(data.sb_r.^2 + data.sb_z.^2);
    br_hat = data.sb_r./bpol;
    bz_hat = data.sb_z./bpol;
    data.krho = kr.*bz_hat - kz.*br_hat;
    data.kpol = kr.*br_hat + kz.*bz_hat;
    data.ktol = kphi;
    rhor = (0.5/pi) * spline(data.rho_bin_center, data.pollen, data.spsi);
    data.wnrho = data.krho*c./w;
    data.mpol = rhor.*data.kpol;
    data.ntol = data.wr.*data.ktol;

    % Calculate the group velocity vgrho, vgpol, vgtol
    data.vgrho = data.vgr_r.*bz_hat - data.vgr_z.*br_hat;
    data.vgpol = data.vgr_r.*br_hat + data.vgr_z.*bz_hat;
    data.vgtol = data.vgr_phi;

    % Calculate damping per species
    if isfield(data, 'salphas')
        % Total damping rate(ki)
        if nargin < 2 && ~isfield(data, 'i_salphal')
            data.i_salphal = ones(size(data.charge));
        elseif nargin >= 2 && ~isfield(data, 'i_salphal')
            data.i_salphal = i_salphal;
        end
        Ind = data.i_salphal ~= 1;
        data.salpha = data.salphal + data.salphac + sum(data.salphas(:, :, Ind), 3); % salphac is likely included in salphal
        Ind1 = data.salpha == 0 | isnan(data.salpha) | isinf(data.salpha);
        % Damping rate(ki) per species
        data.salphae = data.salpha - data.salphac - sum(data.salphas, 3);
        data.salphas = cat(3, data.salphae, data.salphas);
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %
        % sum(delpwr_s) and sum(delpwr_c) may be                  %
        % O(1%) smaller than the GENRAY results                    %
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %
        % Collisional power deposition
        data.delpwr_c = data.deposite.*abs(data.salphac)./data.salpha;
        data.delpwr_c(Ind1) = 0;
        % Power deposition per species
        data.delpwr_s = data.salphas.*abs(data.deposite)./data.salpha;
        for k = 1:size(data.delpwr_s, 3)
            temp = data.delpwr_s(:, :, k);
            temp(Ind1) = 0;
            data.delpwr_s(:, :, k) = temp;
        end
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %
        % sum(delpwrcs) and sum(delpwrcs_c) may be               %
        % O(1%) smaller than the GENRAY results                    %
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %
        % Cross-section power deposition
        data.delpwrcs = zeros([size(data.eqdsk_psi.'), data.nbulk]);
        data.delpwrcs_v = zeros([size(data.eqdsk_psi.'), data.nbulk]);
        data.delpwrcs_c = zeros(size(data.eqdsk_psi.'));
        data.delpwrcs_cv = zeros(size(data.eqdsk_psi.'));
        for j = 1:data.nbulk
            [data.delpwrcs(:, :, j), data.delpwrcs_v(:, :, j)] =...
                GENRAY.ray2cs(data.eqdsk_r.', data.eqdsk_z, data.wr, data.wz, data.delpwr_s(:, :, j), data.nrayelt);
        end
        [data.delpwrcs_c, data.delpwrcs_cv] =...
            GENRAY.ray2cs(data.eqdsk_r.', data.eqdsk_z, data.wr, data.wz, data.delpwr_c, data.nrayelt);
    end

    % Classify rays with Npar and launch position
    data = GENRAY.rayclass(data);
end