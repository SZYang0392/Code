function s = SigmaDT(T)
% Calculate reaction rate (m^3/s) of DT fusion
% Ref :  P.T. Bonoli and M. Porkolab 1987 Nucl. Fusion 27 1341
    TkeV = T/1e3;
    s = 3.68e-18*exp(-20./(TkeV.^(1/3)))./(TkeV.^(2/3));
end