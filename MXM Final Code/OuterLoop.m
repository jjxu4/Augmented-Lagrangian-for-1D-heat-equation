function q_optimal = OuterLoop(q0, params, source, endpoints)

    % Basic parameters
    rho = params.rho_initial;
    mu = params.mu_initial;
    gamma = params.gamma;
    tau = params.tau;
    R = params.R;
    epsilon = params.epsilon;

    u_d   = params.u_d;
    p_d   = params.p_d;
    u_max = params.u_max;
    p_max = params.p_max;

    k = 1; % steps
    n =  1; % successful steps

    q = q0; % initial guess

    while R > epsilon
        % Step 1
        [q_opt, u_opt, p_opt, fval, exitflag, output] = InnerLoop(q, params, source, endpoints);
    
        % TODO: display some diagnostics
    
        % Step 2
        mu_1_potential = max(rho * (u_opt - u_max) + mu, 0);
        mu_2_potential = max(rho * (-u_opt) + mu, 0);
        mu_3_potential = max(rho * (p_opt - p_max) + mu, 0);
        mu_4_potential = max(rho * (-u_opt) + mu, 0);
    
        % Step 3 (INCOMPLETE, need to add second part with the |int( ???)|)
        error_1 = max(abs(u_opt - u_d));
        error_2 = max(abs(-u_opt));
        error_3 = max(abs(p_opt - p_d));
        error_4 = max(abs(-p_opt));
    
        R_next = max([error_1, error_2, error_3, error_4]);
    
        if R_next <= tau * R
            % Successful step
            mu{1} = mu_1_potential;
            mu{2} = mu_2_potential;
            mu{3} = mu_3_potential;
            mu{4} = mu_4_potential;
            n = n + 1;
            disp("Succesful step")
    
        else
            % Unsuccessful step
            rho = rho * gamma;
            disp("Unsuccessful step. Increasing penalty parameter")
        end 

        R = R_next;
        k = k + 1;
    end

    q_optimal = q_opt;
end