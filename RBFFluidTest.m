function X = RBFFluidTest(domain, rho, mu, cell, vessel)
    nx = domain.nx;
    ny = domain.ny;
    nz = domain.nz;
    h = (domain.bounds(:, 2) - domain.bounds(:, 1))' ./ [nx ny nz];

    volume = prod(h);

    Ue = zeros(nx * ny * nz, 3);
    Ff = zeros(nx * ny * nz, 3);
    bump = @(x) sin(pi * x).^4;

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
    [solver, L, handles] = NavierStokesSolver(domain, rho, mu, k);

    xs = linspace(0.5/domain.nx, 1 - 0.5/domain.nx, domain.nx);
    ys = linspace(0, 1 - 1/domain.ny, domain.ny);
    zs = linspace(0.5/domain.nz, 1 - 0.5/domain.nz, domain.nz);
    [XS, YS, ZS] = ndgrid(xs, ys, zs);
    Uy = 0; %bump(XS(:)) .* bump(ZS(:));
    Ue(:, 2) = Uy;
    Fy = 1000; %-mu * (L * Ue(:, 2));
    Ff(:, 2) = Fy;

    % Initial configuration is the reference configuration
    X = cell.shape(cell.data_sites);
    Y = vessel.shape(vessel.data_sites);
    X0 = X;

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

        disp(sprintf('%d: %10g %10g %10g\n', s, max(abs(X - X0))));
        if mod(s, 10) == 0
            %U = cell.id * Ulc;
            % Plot
            subplot(1, 2, 1);
            scatter3(X(:, 1), X(:, 2), X(:, 3), 50, r, 'filled');
            hold on;
            ind = Xlv(:, 1) >= 2^-9;
            scatter3(Xlv(ind, 1), mod(Xlv(ind, 2), 2^-8), Xlv(ind, 3), 25, 'filled');
            %quiver3(Xlc(:, 1), Xlc(:, 2), Xlc(:, 3), 1e1 / 64 * U(:, 1), 1e1 / 64 * U(:, 2), 1e1 / 64 * U(:, 3), 'AutoScale', 'off');
            hold off;
            axis equal;
            xlim(domain.bounds(1, :)); ylim(domain.bounds(2, :)), zlim(domain.bounds(3, :));
            title(sprintf('t = %f', s * k));
            view([-90 0]);

            subplot(1, 2, 2);
            scatter3(X(:, 1), X(:, 2), X(:, 3), 50, r, 'filled');
            hold on;
            scatter3(Xlv(:, 1), mod(Xlv(:, 2), 2^-8), Xlv(:, 3), 25, 'filled');
            %quiver3(Xlc(:, 1), Xlc(:, 2), Xlc(:, 3), 1e1 / 64 * U(:, 1), 1e1 / 64 * U(:, 2), 1e1 / 64 * U(:, 3), 'AutoScale', 'off');
            hold off;
            axis equal;
            xlim(domain.bounds(1, :)); ylim(domain.bounds(2, :)), zlim(domain.bounds(3, :));
            title(sprintf('t = %f', s * k));
            view([0 0]);
            
            drawnow;
        end
    end
end
