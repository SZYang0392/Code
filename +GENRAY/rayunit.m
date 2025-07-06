function C = rayunit(Char)
    switch Char
        case 'ws'
            C = '$s / m$';
        case 'seikon'
            C = '$det(D)$';
        case 'spsi'
            C = '$\rho$';
        case 'wr'
            C = '$R / m$';
        case 'wphi'
            C = '$\zeta / rad$';
        case 'wz'
            C = '$Z / m$';
        case 'w_theta_pol'
            C = '$\theta$ / degree';
        case 'wnpar'
            C = '$N_{\parallel}$';
        case 'wnper'
            C = '$N_{\perp}$';
        case 'delpwr'
            C = '$P / W$';
        case 'sdpwr'
            C = '$k_{i, i}$(all) / $m^{-1}$';
        case 'salphas'
            C = '$k_{i, i} / m^{-1}$';
        case 'wdnpar'
            C = '$dN_{\parallel}$';
        case 'cwexde'
            C = '$E_x / V\cdot m^{-1}$';
        case 'cweyde'
            C = '$E_y / V\cdot m^{-1}$';
        case 'cwezde'
            C = '$E_z / V\cdot m^{-1}$';
        case 'fluxn'
            C=  '$S / W\cdot m^{-2}$';
        case 'sbtot'
            C = '$B / T$';
        case 'sene'
            C = '$n_e / cm^{-3}$';
        case 'ste'
            C = '$T_e / eV$';
        case 'szeff'
            C = '$Z_{eff}$';
        case 'salphac'
            C = '$k_{i, c} / m^{-1}$';
        case 'salphal'
            C = '$k_i / m^{-1}$';
        case 'sb_r'
            C = '$B_R / T$';
        case 'sb_z'
            C = '$B_z / T$';
        case 'sb_phi'
            C = '$B_t / T$';
        case 'wn_r'
            C = '$N_r$';
        case 'wn_z'
            C = '$N_z$';
        case 'wn_phi'
            C = '$N_{\zeta}$';
        case 'vgr_r'
            C = '$v_{gr} / m\cdot s^{-1}$';
        case 'vgr_z'
            C = '$v_{gz} / m\cdot s^{-1}$';
        case 'vgr_phi'
            C = '$v_{g\zeta} / m\cdot s^{-1}$';
        case 'flux_z'
            C = '$S_z / W\cdot m^{-2}$';
        case 'flux_r'
            C = '$S_r / W\cdot m^{-2}$';
        case 'flux_phi'
            C = '$S_{\zeta} / W\cdot m^{-2}$';
        case 'w_eff_nc'
            C = '$J/w / (A\cdot m^{-2}) / (W\cdot m^{-3})$';
        case 'cweps11'
            C = '$\varepsilon_{xx}$';
        case 'cweps12'
            C = '$\varepsilon_{xy}$';
        case 'cweps13'
            C = '$\varepsilon_{xz}$';
        case 'cweps21'
            C = '$\varepsilon_{yx}$';
        case 'cweps22'
            C = '$\varepsilon_{yy}$';
        case 'cweps23'
            C = '$\varepsilon_{yz}$';
        case 'cweps31'
            C = '$\varepsilon_{zx}$';
        case 'cweps32'
            C = '$\varepsilon_{zy}$';
        case 'cweps33'
            C = '$\varepsilon_{zz}$';
        case 'deposite'
            C = '$P / W$';
        case 'salpha'
            C = '$k_{i, total} / m^{-1}$';
        case 'salphae'
            C = '$k_{i, e} / m^{-1}$';
        case 'delpwr_c'
            C = '$P_{coll} / W$';
        case 'delpwr_s'
            C = '$P_{s} / W$';
        otherwise
            C = 'a.u.';
    end
end