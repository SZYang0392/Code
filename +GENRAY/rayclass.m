function genrayout = rayclass(genrayout)
% Classify rays with Npar and launch position

    % -------- Npar -------- %
    Npar = genrayout.wnpar(1, :);
    % Process only finite values (exclude NaN and Inf)
    Npar_finite = Npar(isfinite(Npar));
    
    if isempty(Npar_finite)
        % All values are NaN or Inf: no valid classification
        genrayout.NNpar = 0;
        genrayout.Npar = [];
        genrayout.NparClass = NaN(size(Npar));
        genrayout.NparNum = [];
    else
        genrayout.NparErr = 1e-5;   % Tolerance as originally set
        genrayout.Npar = uniquetol(Npar_finite, genrayout.NparErr);
        genrayout.NNpar = numel(genrayout.Npar);
        % Classify original Npar (including NaN); NaN yields NaN in output
        genrayout.NparClass = discretize(Npar, [genrayout.Npar, Inf]);
        % Count per class; histcounts automatically ignores NaN
        genrayout.NparNum = histcounts(genrayout.NparClass, genrayout.NNpar);
    end

    % -------- Theta -------- %
    Theta = genrayout.w_theta_pol(1, :);
    Theta_finite = Theta(isfinite(Theta));
    
    if isempty(Theta_finite)
        genrayout.Ntheta = 0;
        genrayout.theta = [];
        genrayout.thetaClass = NaN(size(Theta));
        genrayout.thetaNum = [];
    else
        genrayout.thetaErr = 1e-5;
        genrayout.theta = uniquetol(Theta_finite, genrayout.thetaErr);
        genrayout.Ntheta = numel(genrayout.theta);
        genrayout.thetaClass = discretize(Theta, [genrayout.theta, Inf]);
        genrayout.thetaNum = histcounts(genrayout.thetaClass, genrayout.Ntheta);
    end

    % Number of launch positions (keeps original logic: equals number of unique angles)
    genrayout.Nlaunch = numel(genrayout.theta);
end