function fisher_rosenbrock_solver()
% FISHER_ROSENBROCK_SOLVER
% Solves the Fisher-KPP reaction-diffusion equation using a 
% semi-implicit Rosenbrock method (1-stage, L-stable) to handle stiffness.
%
% Equation: du/dt = d^2u/dx^2 + u(1-u) on [0,1]x[0,1]
% BCs: Homogeneous Neumann (Zero flux)
% IC:  u(x,0) = cos^2(pi*x)

    % --- PARAMETERS ---
    m = 50;      % Spatial grid subdivisions
    n = 200;     % Temporal steps
    N_ref = 10000; % Fine-grid reference time steps
    
    fprintf('Running simulation with m=%d spatial steps and n=%d time steps...\n', m, n);

    % --- REFERENCE & APPROXIMATION SOLUTIONS ---
    Solu_ex   = solve_method(m, N_ref);
    uex       = Solu_ex(:, end); % Reference solution at T = 1

    Solu_tau  = solve_method(m, n);
    u_tau     = Solu_tau(:, end);
    
    Solu_tau2 = solve_method(m, 2*n);
    u_tau2    = Solu_tau2(:, end);

    % --- VISUALIZATION ---
    plot_2d_comparison(Solu_tau, Solu_tau2, uex, m, n);
    plot_3d_surface(n, m, Solu_tau);

    % --- ERROR & EMPIRICAL CONVERGENCE ANALYSIS ---
    e_tau  = max(abs(u_tau - uex));
    e_tau2 = max(abs(u_tau2 - uex));
    p_emp  = log(e_tau / e_tau2) / log(2);
    
    fprintf('------------------------------------------------\n');
    fprintf('Error (tau)        : %.4e\n', e_tau);
    fprintf('Error (tau/2)      : %.4e\n', e_tau2);
    fprintf('Empirical Order (p): %.4f\n', p_emp);
    fprintf('------------------------------------------------\n');

    export_results(e_tau, e_tau2, p_emp, m, n);
end

% ---------------------------------------------------------
% --- CORE SOLVER FUNCTION (100% VECTORIZED & SPARSE) ---
% ---------------------------------------------------------
function Solu = solve_method(m, n)
    h   = 1 / m;
    x   = linspace(0, 1, m + 1)';
    tau = 1 / n;
    
    % Initial Condition
    W = (cos(pi .* x)).^2;
    
    % Sparse Laplacian Operator with Homogeneous Neumann BCs
    e = ones(m + 1, 1);
    D = spdiags([e, -2*e, e], -1:1, m + 1, m + 1);
    D(1, 2)     = 2;
    D(m + 1, m) = 2;
    Dh = (1 / h^2) * D;
    
    % Dense pre-allocation for solution history (faster column writes)
    Solu = zeros(m + 1, n + 1);
    Solu(:, 1) = W;
    
    % Sparse Identity Matrix (preserves O(m) tridiagonal solver speed)
    I = speye(m + 1);
    
    for j = 1:n
        % Reaction term G(u) = u(1-u) and vectorized sparse Jacobian J_G = diag(1 - 2u)
        G  = W .* (1 - W);
        JG = spdiags(1 - 2 .* W, 0, m + 1, m + 1);
        
        % Linearly Implicit Rosenbrock Step: (I - (tau/2)*(Dh + JG)) * K = tau * (Dh*W + G)
        A = I - (tau / 2) * (Dh + JG);
        B = tau * (Dh * W + G);
        
        K = A \ B;
        W = W + K;
        Solu(:, j + 1) = W;
    end
end

% ---------------------------------------------------------
% --- HELPER FUNCTIONS ---
% ---------------------------------------------------------
function plot_2d_comparison(Solu_tau, Solu_tau2, uex, m, n)
    x = linspace(0, 1, m + 1);
    figure('Name', '2D Solutions', 'Color', 'w');
    
    subplot(2, 1, 1);
    plot(x, Solu_tau(:, 1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Initial Condition u(x,0) = cos^2(\pi x)');
    hold on;
    plot(x, uex, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Steady Evolution at t = 1');
    title('Initial Profile (t = 0) vs. Evolved Profile (t = 1)');
    legend('Location', 'best'); xlabel('x'); ylabel('u(x,t)'); grid on;
    
    subplot(2, 1, 2);
    plot(x, Solu_tau(:, end), 'r.-', 'LineWidth', 1.2, 'DisplayName', sprintf('Step \\tau (n=%d)', n));
    hold on;
    plot(x, Solu_tau2(:, end), 'b--', 'LineWidth', 1.2, 'DisplayName', sprintf('Step \\tau/2 (2n=%d)', 2*n));
    plot(x, uex, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Reference (N=10,000)');
    title(sprintf('Convergence Comparison at Final Time t = 1 (m = %d)', m));
    legend('Location', 'best'); xlabel('x'); ylabel('u(x,1)'); grid on;
end

function plot_3d_surface(n, m, Solu_tau)
    t = linspace(0, 1, n + 1);
    x = linspace(0, 1, m + 1);
    [T, X] = meshgrid(t, x);
    
    figure('Name', '3D Evolution', 'Color', 'w');
    surf(T, X, Solu_tau);
    xlabel('Time (t)'); ylabel('Space (x)'); zlabel('u(x,t)');
    title('Space-Time Evolution of the Fisher-KPP Equation');
    shading interp; colormap jet; colorbar;
    view(-45, 30);
end

function export_results(e_tau, e_tau2, p_emp, m, n)
    fid = fopen('simulation_results.txt', 'w');
    fprintf(fid, 'Fisher-KPP Equation Simulation Results\n');
    fprintf(fid, 'Spatial Grid m = %d | Time Steps n = %d\n', m, n);
    fprintf(fid, 'Max Error (tau)   : %.4e\n', e_tau);
    fprintf(fid, 'Max Error (tau/2) : %.4e\n', e_tau2);
    fprintf(fid, 'Empirical Order p : %.4f\n', p_emp);
    fclose(fid);
end
