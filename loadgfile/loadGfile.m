function G = loadGfile(filename)
% Read G EQDSK file. All variables in a G EQDSK file are in SI units.
    if ~exist(filename, "file")
        G = [];
        fprintf(['Error : file', filename, ' do not exist!\n']);
        return;
    end

    fid = fopen(filename, "r");
    G.comment = fscanf(fid, '%c', 48);
    eFormat = '%e';
    eNum = 5;
    % Scan scalar parameters line by line
    G.Cutline1 = '---------------------------- Scalar Parameters ----------------------------';
    [G.i3, G.nr, G.nz] = fscanf2each(fid, '%d', 3);
    [G.Rboxlen, G.Zboxlen, G.R0, G.Rmin, G.Z0] = fscanf2each(fid, eFormat, eNum);
    [G.Raxis, G.Zaxis, G.Psi_axis, G.Psi_bound, G.B0] = fscanf2each(fid, eFormat, eNum);
    [G.current, G.Psi_axis, ~, G.Raxis, ~] = fscanf2each(fid, eFormat, eNum);
    [G.Zaxis, ~, G.Psi_bound, ~, ~] = fscanf2each(fid, eFormat, eNum);
    % Scan parameters along the redial profile
    G.Cutline2 = '---------------------------- Profile Parameters ----------------------------';
    G.Psi = linspace(G.Psi_axis, G.Psi_bound, G.nr);
    G.Psit = [];
    % Gfile.r = sqrt((Gfile.Psi_sam - Gfile.Psi_axis)/(Gfile.Psi_bound - Gfile.Psi_axis));
    G.RBt = fscanf(fid, eFormat, G.nr).';
    G.p0 = fscanf(fid, eFormat, G.nr).';
    G.dRBt2dPsi = fscanf(fid, eFormat, G.nr).';
    G.dp0dPsi = fscanf(fid, eFormat, G.nr).';
    G.q = [];
    G.Psirz = fscanf(fid, eFormat, [G.nr, G.nz]).';
    G.q = fscanf(fid, eFormat, G.nr).';
    G.Psit = Psip2t(G.Psi, G.q);%!!!!!!!!!!!!!!!!!!!!!!!!Warning in case of the reverse magnetic shear!
    G.rho = sqrt((G.Psit - min(G.Psit))/(max(G.Psit) - min(G.Psit)));
    % Scan the poloidal flux (nr, nz, Rboxlen, Zboxlen, Rmin, Z0)
    G.Cutline3 = '---------------------------- Cross Section Parameters ----------------------------';
    G.R = G.Rmin + (G.Rboxlen/(G.nr-1))*(1:G.nr);
    G.Z = G.Z0 - 0.5*G.Zboxlen + (G.Zboxlen/(G.nz-1))*(1:G.nz);
    G.Z = G.Z.';
    G.inLCFS = false(G.nz, G.nr);
    dR = G.Rboxlen/(G.nr-1);   %From here the values to be calculated.
    Rmid = (G.R(1:end-1) + G.R(2:end))/2;
    dZ = G.Zboxlen/(G.nz-1);
    Zmid = (G.Z(1:end-1) + G.Z(2:end))/2;
    BR_mid = -(G.Psirz(2:end, :) - G.Psirz(1:end-1, :))/dZ./G.R;
    G.BRrz = interp2(G.R, Zmid, BR_mid, G.R, G.Z, 'spline', nan);
    Bz_mid = (G.Psirz(:, 2:end) - G.Psirz(:, 1:end-1))/dR./Rmid;
    G.Bzrz = interp2(Rmid, G.Z, Bz_mid, G.R, G.Z, 'spline', nan);
    G.Bprz = sqrt(G.BRrz.^2 + G.Bzrz.^2);
    G.Btrz = nan(G.nr, G.nz);
    G.Brz = nan(G.nr, G.nz);
    G.Psitrz = nan(G.nr, G.nz);
    G.RBtrz = nan(G.nr, G.nz);
    G.p0rz = nan(G.nr, G.nz);
    G.qrz = nan(G.nr, G.nz);
    G.rhorz = nan(G.nr, G.nz);
    % Scan the plasma boundary and the limiter
    G.Cutline4 = '---------------------------- Boundary Parameters ----------------------------';
    G.nbound = fscanf(fid, '%d', 1).';
    G.nlimiter = fscanf(fid, '%d', 1).';
    G.Rbound = fscanf(fid, eFormat, [2, G.nbound]);
    G.Zbound = G.Rbound(2, :);
    G.Rbound = G.Rbound(1, :);
    G.Rlimiter = fscanf(fid, eFormat, [2, G.nlimiter]);
    G.Zlimiter = G.Rlimiter(2, :);
    G.Rlimiter = G.Rlimiter(1, :);
    % Calculate Indexes of points in the LCFS
    R1 = repmat(G.R, G.nz, 1);
    Z1 = repmat(G.Z, 1, G.nr);
    G.inLCFS = inpolygon(R1, Z1, G.Rbound, G.Zbound);
    % Calculate values on the cross section
    G.RBtrz = Prof2CS(G.Psi, G.RBt, G.Psirz, G.inLCFS);
    G.Btrz = G.RBtrz./G.R;
    G.Brz = sqrt(G.Bprz.^2 + G.Btrz.^2);
    G.Psitrz = Prof2CS(G.Psi, G.Psit, G.Psirz, G.inLCFS);
    G.p0rz = Prof2CS(G.Psi, G.p0, G.Psirz, G.inLCFS);
    G.qrz = Prof2CS(G.Psi, G.q, G.Psirz, G.inLCFS);
    Ind1 = G.inLCFS & G.qrz < 0;
    G.Bprz(Ind1) = -G.Bprz(Ind1);
    G.rhorz = Prof2CS(G.Psi, G.rho, G.Psirz, G.inLCFS);
    % Close file
    fclose(fid);
end

function varargout = fscanf2each(fid, format, num)
    if nargout == 1
        varargout{1} = fscanf(fid, format, num).';
    else
        varargout = cell(1, nargout).';
        output1 = fscanf(fid, format, num);
        for j = 1:length(output1)
            varargout{j} = output1(j);
        end
    end
end