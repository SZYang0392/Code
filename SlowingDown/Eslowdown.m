function E = Eslowdown(m, v0, vc)
% Calculate energy of isotropic slowing-down distribution
    e = 1.602176565e-19;
    E = sqrt(3)*pi + 6*sqrt(3)*acot(sqrt(3)*vc./(2*v0 - vc)) - 3*log((v0 + vc).^2./(v0.^2 - v0.*vc + vc.^2));
    E = -E.*vc.^2./18 + 0.5*v0.^2;
    E = 2*pi*m*E;
    A0 = (3/4/pi)./log(1 + (v0./vc).^3);
    E = E/e.*A0;
end