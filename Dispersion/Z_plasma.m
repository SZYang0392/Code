function Z = Z_plasma(z)
%--------------Plasma dispersion function--------------%
    % Plasma dispersion function Z(zeta).
    % z is a matrix.
    % At least 6 correct digits when 1e-6 <= |z| <= 1e6.
    C = zeros(size(z));
    indp = imag(z) >= 0;
    C(indp) = cef(z(indp), 50);%N controls the effective digits.
    C(~indp) = 2*exp(-z(~indp).^2) - cef(-z(~indp), 50);
    Z = 1i.*1.772453850905516.*C; % sqrt(pi) = 1.772453850905516
end