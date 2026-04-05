function avg_vals = Contour_avg(R, Z, rho, F, rho0_list)
    % Calculates (\int F dl) / (\int dl) along contours where rho = rho0
    % R: 1xn array
    % Z: mx1 array
    % rho, F: mxn arrays
    % rho0_list: vector of rho0 values you want to evaluate
    
    % Initialize output array
    avg_vals = zeros(length(rho0_list), 1);
    
    for k = 1:length(rho0_list)
        rho0 = rho0_list(k);
        
        % Use contourc to get the contour vertices for rho = rho0
        % C contains [level; num_points] followed by [R; Z] coordinates
        C = contourc(R, Z, rho, [rho0 rho0]);
        
        % Accumulate the integrals over all valid loops for this rho0
        total_F_dl = 0;
        total_dl = 0;
        
        idx = 1;
        while idx < size(C, 2)
            % Extract the number of points for this specific contour segment
            num_pts = C(2, idx);
            
            % Extract R and Z coordinates for the contour
            R_c = C(1, idx+1 : idx+num_pts);
            Z_c = C(2, idx+1 : idx+num_pts);
            
            % Advance the index to the next contour segment
            idx = idx + num_pts + 1;
            
            % Calculate the line segment lengths (dl)
            dR = diff(R_c);
            dZ = diff(Z_c);
            dl = sqrt(dR.^2 + dZ.^2);
            
            % Interpolate the 'F' field onto the contour coordinates
            F_c = interp2(R, Z, F, R_c, Z_c, 'spline');
            
            % Average the interpolated F values at the midpoints of each segment
            F_seg = (F_c(1:end-1) + F_c(2:end)) / 2;
            
            % Perform the numerical integration for this loop
            total_F_dl = total_F_dl + sum(F_seg .* dl);
            total_dl = total_dl + sum(dl);
        end
        
        % Calculate the final ratio for this rho0
        if total_dl > 0
            avg_vals(k) = total_F_dl / total_dl;
        else
            avg_vals(k) = NaN; % Return NaN if the contour doesn't exist
        end
    end
end