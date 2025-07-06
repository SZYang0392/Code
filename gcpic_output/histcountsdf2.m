function [fxy, xsam, ysam] = histcountsdf2(x, y, w, Ndist, xsam0, ysam0)
    % xsam0, ysam0 are optional input arguments.
    % do y
    if nargin > 4
        ysam = ysam0;
        dy = (ysam(end) - ysam(1))/(numel(ysam) - 1);
    else
        miny = min(y);
        maxy = max(y);
        dy = (maxy - miny)/Ndist(end);
        miny = miny - dy/2;
        maxy = maxy + dy/2;
        ysam = linspace(miny, maxy, Ndist(end) + 1);
    end
    Indy = y < (ysam(2:end).') & y >= ysam(1:end-1).';
    % do x
    if nargin > 4
        xsam = xsam0;
        dx = (xsam(end) - xsam(1))/(numel(xsam) - 1);
        minx = min(xsam);
    else
        minx = min(x);
        maxx = max(x);
        dx = (maxx - minx)/Ndist(1);
        minx = minx - dx/2;
        maxx = maxx + dx/2;
        xsam = linspace(minx, maxx, Ndist(1) + 1);
    end
    Indx0 = floor((x - minx)/dx) + 1;
    % Combine x, y
    fxy = zeros(Ndist(1), Ndist(end));
    for k = 1:Ndist(1)
        fxy(k, :) = sum(w(Indx0 == k).*Indy(:, Indx0 == k), 2).';
    end
    fxy = fxy/dx/dy/numel(x);
end
