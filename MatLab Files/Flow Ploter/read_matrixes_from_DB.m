function [ni, nj, num_points_on_airfoil, gamma, rho_inf, Mach_inf, p_inf, a_inf, u_inf, alpha_deg, NACA, x, y, rho, u, v, e, p, a, M, p0, CL, CD] = read_matrixes_from_DB(ID, db_table)    
    index = find(db_table.ID == ID);
    x_array = typecast(db_table.x_2Dmat{index},'double');
    y_array = typecast(db_table.y_2Dmat{index},'double');
    rhos    = typecast(db_table.rho_2Dmat{index},'double');
    us      = typecast(db_table.u_2Dmat{index},'double');
    vs      = typecast(db_table.v_2Dmat{index},'double');
    es      = typecast(db_table.e_2Dmat{index},'double');
    
    ni = db_table.ni(index);
    nj = db_table.nj(index);
    num_points_on_airfoil = db_table.num_points_on_airfoil(index);
    gamma     = db_table.Gamma(index);
    rho_inf   = db_table.density(index);
    Mach_inf  = db_table.Mach_inf(index);
    p_inf     = db_table.environment_pressure(index);
    a_inf     = sqrt(gamma * p_inf / rho_inf);
    u_inf     = Mach_inf * a_inf;
    alpha_deg = db_table.angle_of_attack_deg(index);
    NACA      = db_table.NACA(index);
    CL        = db_table.CL(index);
    CD        = db_table.CD(index);
    
    for i = 1:nj
        rho(i,:) = rhos(1+(i-1)*ni:ni+(i-1)*ni);
    end
    rho;
    
    for i = 1:nj
        u(i,:) = us(1+(i-1)*ni:ni+(i-1)*ni);
    end
    u;
    
    for i = 1:nj
        v(i,:) = vs(1+(i-1)*ni:ni+(i-1)*ni);
    end
    v;
    
    for i = 1:nj
        e(i,:) = es(1+(i-1)*ni:ni+(i-1)*ni);
    end
    e;
    
    for i = 1:nj
        x(i,:) = x_array(1+(i-1)*ni:ni+(i-1)*ni);
    end
    x;
    
    for i = 1:nj
        y(i,:) = y_array(1+(i-1)*ni:ni+(i-1)*ni);
    end
    y;
    
    p  = (gamma - 1).*(e - 0.5.*rho.*(u.^2 + v.^2));
    a = sqrt(gamma.*p./rho);
    M = sqrt(u.^2 + v.^2)./a;
    p0 = (p) .* (1 + (gamma-1) / 2 * M.^2).^(gamma / (gamma-1));
end

