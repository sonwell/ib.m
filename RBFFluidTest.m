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
    [solver, L, handles] = UnsteadyStokesSolver(domain, rho, mu, k);

    xs = linspace(0.5/domain.nx, 1 - 0.5/domain.nx, domain.nx);
    ys = linspace(0, 1 - 1/domain.ny, domain.ny);
    zs = linspace(0.5/domain.nz, 1 - 0.5/domain.nz, domain.nz);
    [XS, YS, ZS] = ndgrid(xs, ys, zs);
    Uy = 0; %bump(XS(:)) .* bump(ZS(:));
    Ue(:, 2) = Uy;
    Fy = 1; %-mu * (L * Ue(:, 2));
    Ff(:, 2) = Fy;

    % Initial configuration is the reference configuration
    X = cell.shape(cell.data_sites);
    Y = vessel.shape(vessel.data_sites);

    for s = 0:150000
        % Backward Euler step
        [Flc, Xlc] = cell.force(X);
        [Flv, Xlv] = vessel.force(Y);
        Fc = spread(Xlc, Flc, cell.ds);
        Fv = spread(Xlv, Flv, vessel.ds);
        [Us, ~] = solver(handles{1}, Ue, Fc + Fv + Ff);
        Ulc = interpolate(X, Us, volume);
        Ulv = interpolate(Y, Us, volume);
        Xi = X + (k/2) * Ulc;
        Yi = Y + (k/2) * Ulv;

        % Crank-Nicolson step
        [Flc, Xlc] = cell.force(Xi);
        [Flv, Xlv] = vessel.force(Yi);
        Fc = spread(Xlc, Flc, cell.ds);
        Fv = spread(Xlv, Flv, vessel.ds);
        [Ue, ~] = solver(handles{2}, Ue, Fc + Fv + Ff);
        Ulc = interpolate(X, Ue, volume);
        Ulv = interpolate(Y, Ue, volume);
        X = X + k * Ulc;
        Y = Y + k * Ulv;

        % Plot
        scatter3(Xlc(:, 1), Xlc(:, 2), Xlc(:, 3), 50, cell.ds, 'filled');
        colorbar;
        axis equal;
        xlim(domain.bounds(1, :)); ylim(domain.bounds(2, :)), zlim(domain.bounds(3, :));
        title(sprintf('t = %f', s * k));
        drawnow;
    end
end
