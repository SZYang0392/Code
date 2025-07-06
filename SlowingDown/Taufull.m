function tau = Taufull(Ps, n, E0)
% Cauculate full slowing-down time
% Ref : Phys. Plasmas 30, 092507 (2023)
    Ec0 = Ec(Ps(n).m, Ps(1), Ps([2:n-1, n+1:numel(Ps)]));
    tause = Tause(Ps, n);
    tau = tause.*log(1 + (E0./Ec0).^(3/2))/3;
end