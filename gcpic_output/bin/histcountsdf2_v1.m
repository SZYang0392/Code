function [fxy, xsam, ysam] = histcountsdf2(x, y, w, Ndist)
    % do x
    minx = min(x);
    maxx = max(x);
    dx = (maxx - minx)/Ndist(1);
    minx = minx - dx/2;
    maxx = maxx + dx/2;
    xsam = linspace(minx, maxx, Ndist(1) + 1);
    Indx = floor((x - minx)/dx);
    % do y
    miny = min(y);
    maxy = max(y);
    dy = (maxy - miny)/Ndist(end);
    miny = miny - dy/2;
    maxy = maxy + dy/2;
    ysam = linspace(miny, maxy, Ndist(end) + 1);
    Indy = floor((y - miny)/dy);
    % Calculate weight
    fxy = zeros(Ndist(1), Ndist(end));
    for j = 1:Ndist(1)
        for k = 1:Ndist(end)
            fxy(j, k) = sum(w(Indx == j & Indy == k));
        end
    end
	fxy = fxy/dx/dy/numel(x);
end
