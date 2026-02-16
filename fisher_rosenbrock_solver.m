function fisher_rosenbrock_solver()
% FISHER_ROSENBROCK_SOLVER
% Solves the Fisher-KPP reaction-diffusion equation using a 
% semi-implicit Rosenbrock method (1-stage) to handle stiffness.
%
% Equation: du/dt = d^2u/dx^2 + u(1-u)
% BCs: Neumann (Zero flux)
% IC: u(x,0) = cos(pi*x)^2

    % --- PARAMETERS ---
    % Grid resolution (Space)
    m = 50; 
    % Time steps
    n = 200; 
    
    fprintf('Running simulation with m=%d spatial steps and n=%d time steps...\n', m, n);

    % --- REFERENCE SOLUTION (High Resolution) ---
    N_ref = 10000;
    Solu_ex = solve_method(m, N_ref);
    uex = Solu_ex(:, end); % Solution at T=1

    % --- APPROXIMATION 1 (Step = tau) ---
    Solu_tau = solve_method(m, n);
    u_tau = Solu_tau(:, end);
    
    % --- APPROXIMATION 2 (Step = tau/2) ---
    Solu_tau2 = solve_method(m, 2*n);
    u_tau2 = Solu_tau2(:, end);

    % --- VISUALIZATION ---
    plot_2d_comparison(Solu_tau, Solu_tau2, uex, m, n);
    plot_3d_surface(n, m, Solu_tau);

    % --- ERROR ANALYSIS (FIXED: Convert sparse to full) ---
    e_tau = full(max(abs(u_tau - uex)));   % <--- FIX AQUÍ
    e_tau2 = full(max(abs(u_tau2 - uex))); % <--- FIX AQUÍ
    
    % Empirical Order of Convergence (p)
    p_emp = log(e_tau / e_tau2) / log(2);
    
    fprintf('------------------------------------------------\n');
    fprintf('Error (tau)   : %.4e\n', e_tau);
    fprintf('Error (tau/2) : %.4e\n', e_tau2);
    fprintf('Empirical Order (p): %.4f\n', p_emp);
    fprintf('------------------------------------------------\n');

    % Save results to file
    export_results(Solu_tau, e_tau, m, n);
end

% ---------------------------------------------------------
% --- CORE SOLVER FUNCTION ---
% ---------------------------------------------------------
function Solu = solve_method(m, n)
    h = 1/m;
    x = linspace(0, 1, m+1)';
    tau = 1/n;
    
    % Initial Condition
    W = (cos(pi .* x)).^2;
    
    % Sparse Matrix Construction (Neumann BCs included)
    e = ones(m+1, 1);
    D = spdiags([e -2*e e], -1:1, m+1, m+1);
    
    % Boundary adjustments (Neumann)
    D(1,2) = D(1,2) + 1;
    D(m+1,m) = D(m+1,m) + 1;
    
    Dh = (1/h^2) * D; % Laplacian Operator
    
    % Pre-allocation
    Solu = sparse(m+1, n+1);
    Solu(:,1) = W;
    
    % Time Integration Loop (Rosenbrock)
    JG = sparse(m+1, m+1);
    I = eye(size(Dh));
    
    for j = 1:n
        % Reaction term G(u) = u(1-u)
        G = W .* (1 - W);
        
        % Jacobian of Reaction Term: J_G = 1 - 2u
        for i = 1:m+1
             JG(i,i) = 1 - 2*W(i);
        end
        
        % Linear System System: (I - tau/2 * J) * K = tau * F(u)
        % J = Dh + JG
        A = I - (tau/2) * (Dh + JG);
        B = tau * (Dh * W + G);
        
        % Solve Linear System
        K = A \ B;
        
        % Update
        W = W + K;
        Solu(:,j+1) = W;
    end
end

% ---------------------------------------------------------
% --- HELPER FUNCTIONS ---
% ---------------------------------------------------------

function plot_2d_comparison(Solu_tau, Solu_tau2, uex, m, n)
    x = linspace(0, 1, m+1);
    figure('Name', '2D Solutions', 'Color', 'w');
    subplot(2,1,1);
    plot(x, Solu_tau(:,1), 'r.-', 'LineWidth', 1.5); hold on;
    plot(x, uex, 'k--', 'LineWidth', 1.5);
    title('Initial Condition (t=0)'); legend('Approx', 'Exact'); grid on;
    
    subplot(2,1,2);
    plot(x, Solu_tau(:,end), 'r.-', 'DisplayName', 'Step \tau'); hold on;
    plot(x, Solu_tau2(:,end), 'b.-', 'DisplayName', 'Step \tau/2');
    plot(x, uex, 'k--', 'DisplayName', 'Reference');
    title(sprintf('Solution at Final Time t=1 (m=%d)', m));
    legend; xlabel('x'); ylabel('u(x,t)'); grid on;
end

function plot_3d_surface(n, m, Solu_tau)
    t = linspace(0, 1, n+1);
    x = linspace(0, 1, m+1);
    [T, X] = meshgrid(t, x);
    
    figure('Name', '3D Evolution', 'Color', 'w');
    surf(T, X, full(Solu_tau)); % Use full() for plotting sparse matrices
    xlabel('Time (t)'); ylabel('Space (x)'); zlabel('u(x,t)');
    title('Time Evolution of Fisher Equation');
    shading interp; colormap jet; colorbar;
    view(-45, 30);
end

function export_results(Solu, error_abs, m, n) 
    % No es necesario guardar todos los puntos si no quieres
    fid = fopen('simulation_results.txt', 'w');
    fprintf(fid, 'Fisher Equation Simulation Results\n');
    fprintf(fid, 'Grid: %dx%d, Max Error: %.4e\n', m, n, full(error_abs));
    fclose(fid);
end