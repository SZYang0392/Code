function w = cef(z,N)
%--------------Complex error function--------------%
    % https://doi.org/10.1137/0731077
    % Computes the function w(z) = exp(-z^2)*erfc(-iz) using a rational
    % series with N terms. It is assumed that Im(z) > 0 or Im(z) = 0.
    % z can be a vector or matrix.
    M = 2*N; M2 = 2*M; k = (-M+1:1:M-1)';
    L = sqrt(N/sqrt(2));
    theta = k*pi/M; t = L*tan(theta/2);
    f = exp(-t.^2).*(L^2+t.^2); f = [0; f];
    a = real(fft(fftshift(f)))/M2;
    a = flipud(a(2:N+1));
    Z = (L + 1i*z)./(L - 1i*z); p = polyval(a,Z);
    w = 2*p./(L - 1i*z).^2 + (1/sqrt(pi))./(L - 1i*z);
end