function L = AugLagrange(q, params, source, endpoints)
% Computes the scalar value of the Augmented Lagrangian
    
    % Helper function
    function f = norm_squared(v)
    % Calculates square norm of vector ( |v|^2 )
        f = dot(v, v);
    end

    % PLACEHOLDER !!!!!
    source.F=q;
    [u, p] = approxPVEsol(params,source,endpoints,2);
    u = u(:).';
    p = p(:).';

    % Basic parameters
    W      = params.W;        % scalar dx*dt
    beta = params.beta;

    u_d   = params.u_d;
    p_d   = params.p_d;
    u_max = params.u_max;
    p_max = params.p_max;

    mu = params.mu;
    rho = params.rho;

    % First half
    u_portion = norm_squared(u - u_d) * W * 0.5;
    p_portion = norm_squared(p - p_d) * W * 0.5;
    q_portion = norm_squared(q) * W * beta * 0.5;

    % Second half
    penalty_1 = (norm_squared(max(mu{1} + rho{1} * (u - u_max), 0)) - norm_squared(mu{1}))* W * 0.5 * (1/rho{1});
    penalty_2 = (norm_squared(max(mu{2} + rho{2} * (-1*u), 0)) - norm_squared(mu{2}))* W * 0.5 * (1/rho{2});
    penalty_3 = (norm_squared(max(mu{3} + rho{3} * (p - p_max), 0)) - norm_squared(mu{3})) * W * 0.5 * (1/rho{3});
    penalty_4 = (norm_squared(max(mu{4} + rho{4} * (-1*p), 0)) - norm_squared(mu{4})) * W * 0.5 * (1/rho{4});

    % return
    L = u_portion + p_portion + q_portion + penalty_1 + penalty_2 + penalty_3 + penalty_4;
end


function f = norm_squared(v)
% Calculates square norm of vector ( |v|^2 )
    f = dot(v, v);
end
