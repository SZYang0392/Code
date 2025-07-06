function [fx, xsam] = histcountsdf(x, w, Ndist)
    minx = min(x);
    maxx = max(x);
    dx = (maxx - minx)/Ndist;
    minx = minx - dx/2;
    maxx = maxx + dx/2;
    xsam = linspace(minx, maxx, Ndist + 1);
    Ind = x < xsam(2:end).' & x >= xsam(1:end-1).';
    fx = sum(w.*Ind, 2).';
	fx = fx/dx/numel(x);
end
