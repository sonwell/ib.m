function X = RBFFluidTest(domain, CP, EP, rho, mu, shape)
    alpha = CP(:, 1);
    beta = CP(:, 2);
    gamma = EP(:, 1);
    delta = EP(:, 2);

    structure = Structure(1, alpha, beta, gamma, delta, shape);

    nx = domain.nx;
    ny = domain.ny;
    nz = domain.nz;
    h = (domain.bounds(:, 2) - domain.bounds(:, 1))' ./ [nx ny nz];

    force = ForceCalculator(structure, 2.5e-3, 2.5e-1, 0);
    sa = SurfaceArea(structure);
    volume = prod(h);

    Ue = zeros(nx * ny * nz, 3);
    Ff = zeros(nx * ny * nz, 3);
    Fy = zeros(nx * ny * nz, 1);
    bump = @(x) sin(pi * x).^4;

    offsets = [0.0, 0.5, 0.5; 0.5, 0.0, 0.5; 0.5, 0.5, 0.0];
    spread = @(x, f, da) ...
        [Transfer('l2e', domain, x, offsets(1, :), f(:, 1) .* da), ...
         Transfer('l2e', domain, x, offsets(2, :), f(:, 2) .* da), ...
         Transfer('l2e', domain, x, offsets(3, :), f(:, 3) .* da)];
    interpolate = @(x, u, dv) ...
        [Transfer('e2l', domain, x, offsets(1, :), u(:, 1) .* dv), ...
         Transfer('e2l', domain, x, offsets(2, :), u(:, 2) .* dv), ...
         Transfer('e2l', domain, x, offsets(3, :), u(:, 3) .* dv)];

    k = 0.015 * exp(mean(log(h)));
    [solver, L, handles] = UnsteadyStokesSolver(domain, rho, mu, k);

    xs = linspace(0.5/domain.nx, 1 - 0.5/domain.nx, domain.nx);
    ys = linspace(0, 1 - 1/domain.ny, domain.ny);
    zs = linspace(0.5/domain.nz, 1 - 0.5/domain.nz, domain.nz);
    [XS, YS, ZS] = ndgrid(xs, ys, zs);
    Uy = bump(XS(:)) .* bump(ZS(:));
    Ue(:, 2) = Uy;
    Fy = -mu * (L * Ue(:, 2));
    Ff(:, 2) = Fy;

    X = shape(alpha, beta);
    for s = 0:400
        [Fl, Xl] = force(X);
        Fc = spread(Xl, Fl, sa);
        [Us, ~] = solver(handles{1}, Ue, Fc + Ff);
        Ul = interpolate(X, Us, volume);
        Xi = X + (k/2) * Ul;

        [Fl, Xl] = force(Xi);
        Fc = spread(Xl, Fl, sa);
        [Ue, ~] = solver(handles{2}, Ue, Fc + Ff);
        Ul = interpolate(X, Ue, volume);
        X = X + k * Ul;

        scatter3(Xl(:, 1), Xl(:, 2), Xl(:, 3), 200, sa, 'filled');
        axis equal;% axis vis3d;
        xlim(domain.bounds(1, :)); ylim(domain.bounds(2, :)), zlim(domain.bounds(3, :));
        title(sprintf('t = %f', s * k));
        drawnow;
    end
end
