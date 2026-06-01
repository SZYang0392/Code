function writeGfile(G, filename, filename_original)
    % Creage a g eqdsk file with original data

    % Open file
    fid = fopen(filename, "w");
    if fid < 0
        fprintf('In writeGfile.m : failed to open target file!\n');
        fprintf(['filename = ', fname, '\n']);
        return;
    end

    % ----------------- Write data ----------------- %
    FortW = 16;
    FortD = 9;
    xdum = 0;
    % 1 : the first line
    fprintf(fid, G.comment);
    fprintf(fid, '   %d %d %d\n', G.i3, G.nr, G.nz);
    % 2 : the scalar parameters
    PrintArray(fid, [G.Rboxlen, G.Zboxlen, G.R0, G.Rmin, G.Z0], FortW, FortD);
    PrintArray(fid, [G.Raxis, G.Zaxis, G.Psi_axis, G.Psi_bound, G.B0], FortW, FortD);
    PrintArray(fid, [G.current, G.Psi_axis, xdum, G.Raxis, xdum], FortW, FortD);
    PrintArray(fid, [G.Zaxis, xdum, G.Psi_bound, xdum, xdum], FortW, FortD);
    % 3 : the profile parameters
    PrintArray(fid, G.RBt, FortW, FortD);
    PrintArray(fid, G.p0, FortW, FortD);
    PrintArray(fid, G.dRBt2dPsi, FortW, FortD);
    PrintArray(fid, G.dp0dPsi, FortW, FortD);
    PrintArray(fid, G.Psirz.', FortW, FortD, false, false);
    PrintArray(fid, G.q, FortW, FortD);
    % 4 : the boundary parameters
    fprintf(fid, '   %d   %d\n', G.nbound, G.nlimiter);
    PrintArray(fid, [G.Rbound; G.Zbound], FortW, FortD);
    PrintArray(fid, [G.Rlimiter; G.Zlimiter], FortW, FortD);

    % ----------------- Optional : copy the remaining of the original file ----------------- %
    if nargin == 3
        fid_original = fopen(filename_original, "r");
        if fid_original < 0
            fprintf('In writeGfile.m : failed to open original file!\n');
            fprintf(['filename_original = ', filename_original, '\n']);
            return;
        end
        ndata = 23 + 5*G.nr + G.nr*G.nz + 2 + 2*(G.nbound + G.nlimiter);
        fscanf(fid_original, '%c', 48);
        fscanf(fid_original, '%e', ndata);
        fscanf(fid_original, '%c', 1);
        remaining = fscanf(fid_original, '%c', Inf);
        fclose(fid_original);
        fprintf(fid, remaining);
    end

    % ----------------- Close file ----------------- %
    fclose(fid);
    if nargin == 3
        fprintf('In writeGfile.m : original file copied successfully!\n');
        fprintf(['filename_original = ', filename_original, '\n']);
    end
end

function PrintArray(fid, array, FortW, FortD, leading_zero, CapitalizeE)
    if nargin < 5
        leading_zero = true;
    end
    if nargin < 6
        CapitalizeE = true;
    end
    iprint = 0;
    for j = 1:numel(array)
        sprint = fortran_ew_d(array(j), FortW, FortD, leading_zero, CapitalizeE);
        fprintf(fid, '%s', sprint);
        iprint = iprint + 1;
        if iprint == 5
            fprintf(fid, '\n');
            iprint = 0;
        end
    end
    if iprint ~= 0
        fprintf(fid, '\n');
    end
end