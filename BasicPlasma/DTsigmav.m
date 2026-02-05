function sigmav = DTsigmav(T_eV)
% Calculate <sigma v> in m^3/s for D-T fusion with
%  L.M. Hively 1977 Nucl. Fusion 17 873 <OR> 
% Freidberg, Jeffrey P. Plasma physics and fusion energy. Cambridge
% university press, 2008. Problem 3.1
% T in eV
    T_keV = T_eV/1e3;
    sigmav = 1e-6 * exp(-21.38.*T_keV.^-0.2935 - 25.20 - 7.101e-2*T_keV ...
        + 1.938e-4*T_keV.^2 + 4.925e-6*T_keV.^3 - 3.984e-8.*T_keV.^4);
end
