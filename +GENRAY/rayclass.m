function genrayout = rayclass(genrayout)
% Classify rays with Npar and launch position
    % Npar
    Npar = genrayout.wnpar(1, :);
    genrayout.NNpar = 0;
    genrayout.NparErr = 1e-5;
    genrayout.Npar = uniquetol(Npar, genrayout.NparErr);
    genrayout.NNpar = numel(genrayout.Npar);
    genrayout.NparClass = discretize(Npar, [genrayout.Npar, Inf]);
    genrayout.NparNum = histcounts(genrayout.NparClass, genrayout.NNpar);
    % Theta
    Theta = genrayout.w_theta_pol(1, :);
    genrayout.Ntheta = 0;
    genrayout.thetaErr = 1e-5;
    genrayout.theta = uniquetol(Theta, genrayout.thetaErr);
    genrayout.Ntheta = numel(genrayout.theta);
    genrayout.thetaClass = discretize(Theta, [genrayout.theta, Inf]);
    genrayout.thetaNum = histcounts(genrayout.thetaClass, genrayout.Ntheta);
    % Launch position
    genrayout.Nlaunch = numel(genrayout.theta);
end