function q_optimal = OuterLoop(q_initial, params, source, endpoints)
% Implements Algorithm 2 in 2024 paper. 

% q_initial is expected to be an (nx + 1) X (nt) grid
% q_n is returned as an (nx+1) X (nt) grid

    % Basic parameters
    rho = params.rho_initial;
    mu = params.mu_initial;
    gamma = params.gamma;
    tau = params.tau;
    R_0 = params.R;
    epsilon = params.epsilon;   % this is the optimality tolerance

    nx = params.nx;
    nt = params.nt;

    W = params.W;

    % history of residuals and penalty parameters (indexed by successful steps)
    R_hist  = [];           
    rho_hist = [];              

    % state constraints. These are expected to be vectors
    u_max = params.u_max;
    p_max = params.p_max;
    u_min = params.u_min;
    p_min = params.p_min;

    % Force everything to be column vectors
    u_max = u_max(:); u_min = u_min(:);
    p_max = p_max(:); p_min = p_min(:);
    for i = 1:4
        mu{i} = mu{i}(:);
    end


    k = 1; % steps
    n =  1; % successful steps

    q_n = q_initial; % initialize guess
    R_n_minus_1 = R_0;
    disp("======================== Starting Outer Loop ========================")
    disp("=====================================================================")
    while R_n_minus_1 > epsilon
        fprintf('\n--------- Outer Iteration %d ---------\n', k);
        fprintf('Current residual R_{k-1} = %.6e\n', R_n_minus_1);


        % Step 1
        params.mu  = mu;
        params.rho = rho;
        [q_k_bar, u_k_bar, p_k_bar, fval, exitflag, output] = InnerLoop(q_n, params, source, endpoints);
        
        % InnerLoop returns all as grids. we need to convert them to
        % vectors
        q_k_bar = q_k_bar(:);
        u_k_bar = u_k_bar(:);
        p_k_bar = p_k_bar(:);

        % Step 2
        mu_1_bar = max(rho{1} * (u_k_bar - u_max) + mu{1}, 0);
        mu_2_bar = max(rho{2} * (u_min - u_k_bar) + mu{2}, 0);
        mu_3_bar = max(rho{3} * (p_k_bar - p_max) + mu{3}, 0);
        mu_4_bar = max(rho{4} * (p_min - p_k_bar) + mu{4}, 0);
    
        % Step 3 
        error_1 = max( max((u_k_bar - u_max), 0) ) + abs( W * sum( mu_1_bar .* (u_max - u_k_bar) ) );
        error_2 = max( max((u_min - u_k_bar), 0) ) + abs( W * sum( mu_2_bar .* (u_k_bar - u_min) ) );
        error_3 = max( max((p_k_bar - p_max), 0) ) + abs( W * sum( mu_3_bar .* (p_max - p_k_bar) ) );
        error_4 = max( max((p_min - p_k_bar), 0) ) + abs( W * sum( mu_4_bar .* (p_k_bar - p_min) ) );
    
        R_k = max([error_1, error_2, error_3, error_4]);
    
        % Step 4
        if R_k <= tau * R_n_minus_1
            % Successful step
        
            % update mu
            mu{1} = mu_1_bar;
            mu{2} = mu_2_bar;
            mu{3} = mu_3_bar;
            mu{4} = mu_4_bar;
        
            % rho remains the same
            % update optimal q. need to convert back into grid, because
            % that is what InnerLoop expects
            q_n = q_k_bar;
            q_n = reshape(q_n, nx+1, nt);
        
            % update R
            R_n_minus_1 = R_k;
        
            % log history for this successful iteration
            R_hist(end+1,1)     = R_k;
            rho_hist(end+1,1:4) = [rho{1}, rho{2}, rho{3}, rho{4}];
        
            % increase succesful step counter
            n = n + 1;
            disp("Succesful step")

    
        else
            % Unsuccessful step
            
            % mu remains unchanged
            % Increase penalty parameter
            rho{1} = rho{1} * gamma;
            rho{2} = rho{2} * gamma;
            rho{3} = rho{3} * gamma;
            rho{4} = rho{4} * gamma;
            disp("Unsuccessful step. Increasing penalty parameter")

            % use a warm guess for q in the next iteration. Need to reshape
            % because that is what InnerLoop expects
            q_n = q_k_bar;
            q_n = reshape(q_n, nx+1, nt);
        end 
        % update counter
        k = k + 1;
    end
    disp('--------- Outer Loop Finished ---------');
    fprintf('Final residual R_{k-1} = %.6e\n', R_n_minus_1);
    q_optimal = q_n;

    % ===================== Plot histories ==========================
    if ~isempty(R_hist)
        succ_iter = 1:numel(R_hist);

        % Plot residual history (log scale is usually informative)
        figure(100);
        semilogy(succ_iter, R_hist, '-o', 'LineWidth', 1.5);
        grid on;
        xlabel('Successful outer iteration');
        ylabel('Residual R_k');
        title('Residual history (successful iterations)');

        % Plot penalty parameter history (log scale)
        figure(101);
        semilogy(succ_iter, rho_hist(:,1), '-o', ...
                 succ_iter, rho_hist(:,2), '-s', ...
                 succ_iter, rho_hist(:,3), '-^', ...
                 succ_iter, rho_hist(:,4), '-d', 'LineWidth', 1.5);
        grid on;
        xlabel('Successful outer iteration');
        ylabel('\rho');
        title('Penalty parameters history (successful iterations)');
        legend({'\rho_1','\rho_2','\rho_3','\rho_4'}, 'Location', 'best');

    end
    % ===============================================================
end

