function varargout = ray2cs(R, Z, rayr, rayz, rayf, nrayelt)
    % Clear values after the end of rays
    Size = size(rayf);
    Size(1) = 1;
    Ind0 = repmat((1:size(rayf, 1)).', Size);
    nrayelt = reshape(nrayelt, 1, []);
    Ind0 = Ind0 < nrayelt;
    % Reshape to 1-D array
    rayr = reshape(rayr(Ind0), [], 1);
    rayz = reshape(rayz(Ind0), [], 1);
    rayf = reshape(rayf(Ind0), [], 1);
    rayf(isnan(rayf) | isinf(rayf)) = 0;
    % Mech indices
    [~, ~, IndZ] = histcounts(rayz, Z);
    IndZ(IndZ >= numel(Z)) = numel(Z) - 1;
    zcoeff = (rayz - Z(IndZ))./(Z(IndZ+1) - Z(IndZ));
    zcoeff = zcoeff(:);
    [~, ~, IndR] = histcounts(rayr, R);
    IndR(IndR >= numel(R)) = numel(R) - 1;
    rcoeff = (rayr - R(IndR).')./(R(IndR+1).' - R(IndR).');
    rcoeff = rcoeff(:);
    % wzr
    wnn = (1 - zcoeff).*(1 - rcoeff);
    wpn = zcoeff.*(1 - rcoeff);
    wnp = (1 - zcoeff).*rcoeff;
    wpp = zcoeff.*rcoeff;
    % Accumulate ray power
    Ind = [IndZ, IndR; IndZ+1, IndR; IndZ, IndR+1; IndZ+1, IndR+1];
    Value = [wnn.*rayf; wpn.*rayf; wnp.*rayf; wpp.*rayf];
    Zmin = min(Z); Zmax = max(Z);
    Rmin = min(R); Rmax = max(R);
    valid = (rayz>=Zmin) & (rayz<=Zmax) & (rayr>=Rmin) &(rayr<=Rmax);
    valid = repmat(valid, 4, 1);
    Value(~valid) = 0;
    Dcs = accumarray(Ind, Value, [numel(Z), numel(R)]);
    % Calculate density
    dZ = (Z(3:end) - Z(1:end-2))/2;
    dZ = [Z(2) - Z(1); dZ; Z(end) - Z(end-1)];
    Rmid = (R(2:end) + R(1:end-1))/2;
    Rmid = [Rmid(1)-R(2)+R(1), Rmid, Rmid(end)+R(end)-R(end-1)];
    dV = pi.*(Rmid(2:end).^2 - Rmid(1:end-1).^2).*dZ;
    Dcs_v = Dcs./dV;

    varargout{1} = Dcs;
    varargout{2} = Dcs_v;
end