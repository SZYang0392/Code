function Npar = Npar_Cardinali2007(w, B, Ps, R0, a, NparA, Theta0, q0)
% Calculate Npar upshift analytically
    %% Physical Constants
    c = 299792458;
    me = 9.1093829e-31;
    e = 1.602176565e-19;
    u = 1.660538921e-27;
    epsilon_0 = 8.854187817e-12;

    %% Basic Parameters
    epsi = a./R0;
    delta0 = c./w./a;
    w2 = w.^2;
    wpi2 = zeros(size(Ps(2).n0));
    for k = 2:numel(Ps)
        Ps(k).wp2 = (Ps(k).q.*e).^2.*Ps(k).n0/epsilon_0./Ps(k).m;
        Ps(k).wc = abs(Ps(k).q).*e.*B./Ps(k).m;
        if k >= 2
            wpi2 = wpi2 + Ps(k).wp2;
        end
    end
    wpi2 = wpi2./w2;
    wpi20 = wpi2(1);
    wpe2 = Ps(1).wp2./w2;
    wpe20 = wpe2(1);
    delta = wpi20/wpe20;
    lambda = w2./((wce.*delta).^2);

    %% Initial Parameters
    Npar0 = NparA.*(1 + epsi.*cos(Theta0));

    %% Solution Parameters
    C3 = 1./(delta0.*delta.*epsi).*q0./(wpi20.*lambda).*NparA.*(1 + epsi.*cos(Theta0));
    C4 = sign(NparA).*sqrt(delta).*epsi.^2./q0.*wpi20.^1.5.*lambda;
end