function T = Tsd(Ebth, Ec0)
    Tsd2Ec = (Ebth./Ec0 - 2*g(sqrt(Ebth./Ec0)))./log(1 + (Ebth./Ec0).^(3/2));
    T = Tsd2Ec.*Ec0;
end

function G = g(x)
    G = (1/sqrt(3)).*atan((2.*x - 1)./sqrt(3)) + (1/6).*log((x.^2 - x + 1)./(x + 1).^2) + sqrt(3*pi)/18;
end