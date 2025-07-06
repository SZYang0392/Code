function [fx, xsam] = histcountsdf(x, w, Ndist)
    minx = min(x);
    maxx = max(x);
    dx = (maxx - minx)/Ndist;
    minx = minx - dx/2;
    maxx = maxx + dx/2;
    xsam = linspace(minx, maxx, Ndist + 1);
    Ind = floor((x - minx)/dx);
    fx = zeros(1, Ndist);
    for k = 1:Ndist
        fx(k) = fx(k) + sum(w(Ind == k));
    end
	fx = fx/dx/numel(x);
end
