function X = RBFFluidTest(domain, rho, mu, cell, vessel)
    nx = domain.nx;
    ny = domain.ny;
    nz = domain.nz;
    h = (domain.bounds(:, 2) - domain.bounds(:, 1))' ./ [nx ny nz];

    volume = prod(h);

    Ue = zeros(nx * ny * nz, 3);
    Ff = zeros(nx * ny * nz, 3);

    offsets = [0.0, 0.5, 0.5; 0.5, 0.0, 0.5; 0.5, 0.5, 0.0];  % MAC offsets
    % Helpers functions for spreading and interpolation
    spread = @(x, f, da) ...
        [Transfer('l2e', domain, x, offsets(1, :), f(:, 1) .* da), ...
         Transfer('l2e', domain, x, offsets(2, :), f(:, 2) .* da), ...
         Transfer('l2e', domain, x, offsets(3, :), f(:, 3) .* da)];
    interpolate = @(x, u, dv) ...
        [Transfer('e2l', domain, x, offsets(1, :), u(:, 1) .* dv), ...
         Transfer('e2l', domain, x, offsets(2, :), u(:, 2) .* dv), ...
         Transfer('e2l', domain, x, offsets(3, :), u(:, 3) .* dv)];

    k = 1e-6; %0.015 * exp(mean(log(h)));
    [solver, ~, handles] = NavierStokesSolver(domain, rho, mu, k);

    Uy = 0;
    Ue(:, 2) = Uy;
    Fy = 8000; %-mu * (L * Ue(:, 2));
    Ff(:, 2) = Fy;

    % Initial configuration is the reference configuration
    X = cell.shape(cell.data_sites);
    Y = vessel.shape(vessel.data_sites);

    ls = domain.bounds(:, 2) - domain.bounds(:, 1);
    lw = domain.bounds(:, 1);
    avg = @(X, w) (w' * X) ./ sum(w);
    mx = @(X) [mod(X(:, 1) - lw(1), ls(1)) + lw(1) ...
               mod(X(:, 2) - lw(2), ls(2)) + lw(2) ...
               mod(X(:, 3) - lw(3), ls(3)) + lw(3)];

    fig = figure;
    for s = 0:150000
        % Backward Euler step
        [Flc, Xlc] = cell.force(X);
        [Flv, Xlv] = vessel.force(Y);
        Fc = spread(Xlc, Flc, cell.ds);
        Fv = spread(Xlv, Flv, vessel.ds);
        [Us, ~] = solver(handles{1}, Ue, Ue, Fc + Fv + Ff);
        Ulc = interpolate(X, Us, volume);
        Ulv = interpolate(Y, Us, volume);
        Xi = X + (k/2) * Ulc;
        Yi = Y + (k/2) * Ulv;

        % Crank-Nicolson step
        [Flc, Xlc] = cell.force(Xi);
        [Flv, Xlv] = vessel.force(Yi);
        Fc = spread(Xlc, Flc, cell.ds);
        Fv = spread(Xlv, Flv, vessel.ds);
        [Ue, ~] = solver(handles{2}, Ue, Us, Fc + Fv + Ff);
        Ulc = interpolate(X, Ue, volume);
        Ulv = interpolate(Y, Ue, volume);
        X = X + k * Ulc;
        Y = Y + k * Ulv;

        fprintf('%d: %12g %12g %12g %12g %12g %12g\n', s, avg(Xlc, cell.ds), avg(Flv, vessel.ds));
        if mod(s, 10) == 0
            [csx, csy, csz] = cell.surf(X);
            [vsx, vsy, vsz] = vessel.surf(Y);
            % Plot
            subplot(1, 2, 1);
            ezsurf(csx, csy, csz, [0, pi, 0, 2*pi]);
            hold on;
            ezsurf(vsx, vsy, vsz, [-pi/2, pi/2, 0, 1]);
            hold off;
            axis equal;
            xlim(domain.bounds(1, :)); ylim(domain.bounds(2, :)), zlim(domain.bounds(3, :));
            title(sprintf('$t = %f$', s * k), 'Interpreter', 'latex');
            view([-90 0]);

            subplot(1, 2, 2);
            ezsurf(csx, csy, csz, [0, pi, 0, 2*pi]);
            hold on;
            ezsurf(vsx, vsy, vsz, [0, 2*pi, 0, 1]);
            hold off;
            axis equal;
            xlim(domain.bounds(1, :)); ylim(domain.bounds(2, :)), zlim(domain.bounds(3, :));
            title(sprintf('$t = %f$', s * k), 'Interpreter', 'latex');
            view([0 0]);
            
            drawnow;

            set(fig, 'Units', 'Inches');
            pos = get(fig, 'Position');
            set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
            saveas(fig, sprintf('simulation%03d.pdf', s/1000));
        end
    end
    close all;
end
