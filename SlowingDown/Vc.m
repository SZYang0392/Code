function v = Vc(Ps, Pi)
% Calculate critical velocity of energetic ions
    if nargin < 2
        Pe = Ps(1);
        Pi = Ps(2:end);
    else
        Pe = Ps;
    end
    e = 1.602176565e-19;

    v = zeros(size(Pe.n0));
    vte = sqrt(2*Pe.T*e/Pe.m);
    for k = 1:numel(Pi)
        v = v + (Pe.m/Pi(k).m).*(Pi(k).n0./Pe.n0).*Pi(k).q^2;
    end
    v = (3*sqrt(pi)/4*v).^(1/3);
    v = v.*vte;
end