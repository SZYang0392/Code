function E = Ec(M, Ps, Pi)
% Calculate critical energy of energetic ions
    if nargin < 3
        Pe = Ps(1);
        Pi = Ps(2:end);
    else
        Pe = Ps;
    end
    e = 1.602176565e-19;
    vc = Vc(Pe, Pi);
    E = 0.5*M*vc.^2/e;
end