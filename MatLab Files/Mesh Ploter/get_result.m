function result = get_result(NACA)
    ni                    = 91;
    nj                    = 38;
    num_points_on_airfoil = 61;
    delta_y               = 0.005;
    XSF                   = 1.1;
    YSF                   = 1.1;
    r                     = 0.001;
    omega                 = 0.1;

    result.metadata   = readtable(sprintf('results\\NACA%d\\ni%d\\nj%d\\num_points_on_airfoil%d\\delta_y%g\\XSF%g\\YSF%g\\r%g\\omega%g\\metadata.txt', NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega));
    result.Ls_valuse  = readmatrix(sprintf('results\\NACA%d\\ni%d\\nj%d\\num_points_on_airfoil%d\\delta_y%g\\XSF%g\\YSF%g\\r%g\\omega%g\\Ls_valuse.txt', NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega));
    result.x_mat      = readmatrix(sprintf('results\\NACA%d\\ni%d\\nj%d\\num_points_on_airfoil%d\\delta_y%g\\XSF%g\\YSF%g\\r%g\\omega%g\\x_mat.txt', NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega));
    result.x_mat_init = readmatrix(sprintf('results\\NACA%d\\ni%d\\nj%d\\num_points_on_airfoil%d\\delta_y%g\\XSF%g\\YSF%g\\r%g\\omega%g\\x_mat_init.txt', NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega));
    result.y_mat  = readmatrix(sprintf('results\\NACA%d\\ni%d\\nj%d\\num_points_on_airfoil%d\\delta_y%g\\XSF%g\\YSF%g\\r%g\\omega%g\\y_mat.txt', NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega));
    result.y_mat_init = readmatrix(sprintf('results\\NACA%d\\ni%d\\nj%d\\num_points_on_airfoil%d\\delta_y%g\\XSF%g\\YSF%g\\r%g\\omega%g\\y_mat_init.txt', NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega));
    
end

