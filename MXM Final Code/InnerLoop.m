function [q_opt, u_opt, p_opt, fval, exitflag, output] = InnerLoop(q0, params)
    % Minimizes the Augmented Lagrangian w.r.t q, where q0 is the initial guess

    % Helper function
    function [L, g] = AugLag_and_grad(q, params)
        % Augmented Lagrangian Value at q
        L = AugLagrange(q, params);
    
        % Gradient at q
        % PLACEHOLDER !!!!!
        g = AugLagrangeGrad(q, params);
    end

    % Settings for fminunc function
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'SpecifyObjectiveGradient', true, 'Display', 'iter');

    % returns value of augmented Lagrangian and gradient of augmented Lagrangian
    fun = @(q) AugLag_and_grad(q, params);

    % optimize
    [q_opt, fval, exitflag, output] = fminunc(fun, q0, options);

    % recover optimal states
    [u_opt, p_opt] = approxPVEsol(q_opt, params);
end

