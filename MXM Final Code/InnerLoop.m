function [q_opt, u_opt, p_opt, fval, exitflag, output] = InnerLoop(q0, params, source, endpoints)
    % Minimizes the Augmented Lagrangian w.r.t q, where q0 is the initial guess

    % Helper function
    function [L, g] = AugLag_and_grad(q, params, source, endpoints)
        % Augmented Lagrangian Value at q
        L = AugLagrange(q, params, source, endpoints);
    
        % Gradient at q
        % PLACEHOLDER !!!!!
        g = AugLagrangeGrad(q, source, params, endpoints, params.u_d, params.p_d, params.u_max, params.p_max);
    end

    % Settings for fminunc function
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'SpecifyObjectiveGradient', true, 'Display', 'iter');

    % returns value of augmented Lagrangian and gradient of augmented Lagrangian
    fun = @(q) AugLag_and_grad(q, params, source, endpoints);

    % optimize
    [q_opt, fval, exitflag, output] = fminunc(fun, q0, options);

    % recover optimal states
    [u_opt, p_opt] = approxPVEsol(params, source, endpoints, 2);
end

