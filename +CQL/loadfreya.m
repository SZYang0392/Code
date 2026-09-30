function data = loadfreya(filename)
%CQL.loadfreya  Load FREYA NBI birth points from a freya_points.txt file.
%
%   data = CQL.loadfreya() reads 'freya_points.txt' in the current
%   directory and returns a struct with fields converted to SI:
%     X, Y, Z, R   -- birth position (m)
%     vx, vy, vz   -- birth velocity components (m/s)
%     index        -- birth point sequence number
%     npoints      -- total number of birth points
%
%   data = CQL.loadfreya(filename) uses the specified file.
%
%   All fields are row vectors (1 × npoints).  Subscripting data.R(j),
%   data.Z(j), data.vx(j), etc. all refer to the same particle j.
%
%   The coordinate system (raw CGS from FREYA / CQL3D frplteq.f, converted):
%     X, Y  -- Cartesian coordinates in the horizontal (equatorial) plane
%     Z     -- vertical coordinate (height above the tokamak midplane)
%     R     -- major radius: R = sqrt(X^2 + Y^2)
%     vx,vy,vz -- Cartesian velocity components
%
%   Example:
%       pts = CQL.loadfreya('CFETR_NB_transport_new/freya_points.txt');
%       scatter3(pts.X, pts.Y, pts.Z, 1, sqrt(pts.vx.^2+pts.vy.^2+pts.vz.^2));
%       xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
%       title('FREYA NBI Birth Points coloured by |v|');

    % Default filename
    if nargin < 1 || isempty(filename)
        filename = 'freya_points.txt';
    end

    % Open and read the text file
    fid = fopen(filename, 'r');
    if fid < 0
        error('CQL:loadfreya:FileNotFound', ...
              'Cannot open file: %s', filename);
    end

    % Skip header line
    header = fgetl(fid);

    % Read all numeric data: pnt x y Z R vx vy vz
    % Format in frplteq.f: i7,1x,7ES12.4E2
    C = textscan(fid, '%f%f%f%f%f%f%f%f');
    fclose(fid);

    if isempty(C{1})
        error('CQL:loadfreya:EmptyFile', ...
              'No birth point data found in: %s', filename);
    end

    % Assign fields as row vectors (still CGS)
    data.index   = C{1}(:).';    % row vector
    data.X       = C{2}(:).';
    data.Y       = C{3}(:).';
    data.Z       = C{4}(:).';
    data.R       = C{5}(:).';
    data.vx      = C{6}(:).';
    data.vy      = C{7}(:).';
    data.vz      = C{8}(:).';
    data.npoints = numel(data.index);

    % ---- CGS -> SI conversion -----------------------------------------
    %   lengths:  cm   -> m      (* 0.01)
    %   velocities: cm/s -> m/s  (* 0.01)
    data.X  = data.X  .* 0.01;
    data.Y  = data.Y  .* 0.01;
    data.Z  = data.Z  .* 0.01;
    data.R  = data.R  .* 0.01;
    data.vx = data.vx .* 0.01;
    data.vy = data.vy .* 0.01;
    data.vz = data.vz .* 0.01;

    % Store metadata
    data.filename = filename;
    data.units    = 'SI (m, m/s)';
    data.header   = strtrim(header);

    clear fid header C;
end