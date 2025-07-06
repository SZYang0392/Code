function data = rayclass(data)
% Classify rays with Npar and launch position
    % Npar
    Npar = data.wnpar(1, :);
    data.NNpar = 0;
    data.NparErr = 1e-5;
    data.Npar = uniquetol(Npar, data.NparErr);
    data.NNpar = numel(data.Npar);
    data.NparClass = discretize(Npar, [data.Npar, Inf]);
    data.NparNum = histcounts(data.NparClass, data.NNpar);
    % Theta
    Theta = data.w_theta_pol(1, :);
    data.Ntheta = 0;
    data.thetaErr = 1e-5;
    data.theta = uniquetol(Theta, data.thetaErr);
    data.Ntheta = numel(data.theta);
    data.thetaClass = discretize(Theta, [data.theta, Inf]);
    data.thetaNum = histcounts(data.thetaClass, data.Ntheta);
    % Launch position
    data.Nlaunch = numel(data.theta);
end