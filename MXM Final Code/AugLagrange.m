function L = AugLagrange(q, params, source, endpoints)
% Computes the scalar value of the Augmented Lagrangian
% q is expected to be an (nx + 1) X (nt) grid
% mu{1}, ..., mu{4} are expected to be column vectors (IMPORTANT)
    
    % Helper function
    function f = norm_squared(v)
    % Calculates square norm of vector ( |v|^2 )
        v = v(:);
        f = dot(v, v);
    end

    % We first usue forward solve to retrieve u(q) and p(q).
    source.F = q;
    [u, p] = approxPVEsol(params, source, endpoints, 2);

    % approxPVEsol returns u and p as grids, so we need to convert them
    % into column vectors
    u_vec = u(:);
    p_vec = p(:);
    q_vec = q(:);

    % Parameters for calculating the augmented Lagrangian
    W = params.W;        % scalar dx * dt
    beta = params.beta;
    
    % Convert everything into column vectors for safety
    % Desired solution. 
    u_d   = params.u_d; u_d = u_d(:);
    p_d   = params.p_d; p_d = p_d(:);

    % State constraints
    u_max = params.u_max; u_max = u_max(:);
    p_max = params.p_max; p_max = p_max(:);

    u_min = params.u_min;
    u_min = u_min(:);
    p_min = params.p_min;
    p_min = p_min(:);

    % Outerloop parameters
    mu = params.mu; % cell containing mu{1}, ..., mu{4}
    for i = 1:4
        mu{i} = mu{i}(:);
    end
    
    rho = params.rho; % cell containing rho{1}, ..., rho{4}

    % First half
    u_portion = norm_squared(u_vec - u_d) * W * 0.5;
    p_portion = norm_squared(p_vec - p_d) * W * 0.5;
    q_portion = norm_squared(q_vec) * W * beta * 0.5;

    % Second half
    penalty_1 = (norm_squared(max(mu{1} + rho{1} * (u_vec - u_max), 0)) - norm_squared(mu{1}))* W * 0.5 * (1/rho{1});
    penalty_2 = (norm_squared(max(mu{2} + rho{2} * (u_min - u_vec), 0)) - norm_squared(mu{2}))* W * 0.5 * (1/rho{2});
    penalty_3 = (norm_squared(max(mu{3} + rho{3} * (p_vec - p_max), 0)) - norm_squared(mu{3})) * W * 0.5 * (1/rho{3});
    penalty_4 = (norm_squared(max(mu{4} + rho{4} * (p_min - p_vec), 0)) - norm_squared(mu{4})) * W * 0.5 * (1/rho{4});

    % return
    L = u_portion + p_portion + q_portion + penalty_1 + penalty_2 + penalty_3 + penalty_4;
end