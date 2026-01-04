function x = iterative_reweighted_L1(A, b, N)
    % Parameters
    epsilon = 0.5;           % Given epsilon value
    tolerance = 1e-4;        % Termination tolerance
    max_iter = 1000;         % Maximum number of iterations to prevent infinite loop

    % Initialize W and x
    W = eye(N^3);            % Identity matrix of size N^3 x N^3
    x_prev = zeros(N^3, 1);  % Initial guess of x

    for iter = 1:max_iter
        % Solve the optimization problem: min ||W * x||_1 s.t. A * x = b
        cvx_begin quiet
            variable x(N^3)
            minimize(norm(W * x, 1))
            subject to
                A * x == b
        cvx_end

        % Check for convergence
        if norm(x - x_prev, 2) < tolerance
            break;
        end

        % Update W for next iteration
        for k = 1:N^3
            W(k, k) = 1 / (abs(x(k)) + epsilon);
        end

        % Update x_prev
        x_prev = x;
    end

    % Display results
    fprintf('Algorithm converged in %d iterations\n', iter);
end
