function [solver, L, handles] = UnsteadyStokesSolver(domain, rho, mu, k)
    n = [domain.nx domain.ny domain.nz];
    h = (domain.bounds(:, 2) - domain.bounds(:, 1))' ./ n;

    multiplex = @(fn, n, h) deal(fn(n(1), h(1)), fn(n(2), h(2)), fn(n(3), h(3)));
    identity = @(n, h) speye(n);
    differential = @(n, h) spdiags(ones(n, 1) * [1 -1 1] / h, [-n+1, 0, 1], n, n);
    laplacian = @(n, h) spdiags(ones(n, 1) * [1 1 -2 1 1] / h^2, [-n+1, -1, 0, 1, n-1], n, n);

    [Gx, Gy, Gz] = multiplex(differential, n, h);
    [Lx, Ly, Lz] = multiplex(laplacian, n, h);
    [Ix, Iy, Iz] = multiplex(identity, n, h);

    I = kron(Iz, kron(Iy, Ix));
    L = kron(Lz, kron(Iy, Ix)) + kron(Iz, kron(Ly, Ix)) + kron(Iz, kron(Iy, Lx));
    Dx = kron(Iz, kron(Iy, Gx));
    Dy = kron(Iz, kron(Gy, Ix));
    Dz = kron(Gz, kron(Iy, Ix));
    nL = -L;
    
    lambda = mu * k / (2 * rho);
    H = I - lambda * L;
    C = ichol(H);
    C2 = ichol(nL);

    function du = step(k, u, f)
        rhs = k / rho * (mu * (L * u) + f);

        du = zeros(size(u));
        [dux, flx, rr, itx, rv] = pcg(H, rhs(:, 1), 1e-8, 100000, C, C');
        [duy, fly, rr, ity, rv] = pcg(H, rhs(:, 2), 1e-8, 100000, C, C');
        [duz, flz, rr, itz, rv] = pcg(H, rhs(:, 3), 1e-8, 100000, C, C');
        du(:, 1) = dux;
        du(:, 2) = duy;
        du(:, 3) = duz;
    end

    function [u, p] = update(k, u, du)
        u_star = u + du;
        div_u_star = Dx * u_star(:, 1) + Dy * u_star(:, 2) + Dz * u_star(:, 3);
        pr = div_u_star - mean(div_u_star);
        [dtphi, fl, rr, it, rv] = pcg(nL, -p, 1e-8, 100000, C2, C2');
        u = u_star + [Dx' * dtphi, Dy' * dtphi, Dz' * dtphi];
        p = H * dtphi / k;
    end

    function [u, p] = be_step(u, f)
        du = step(k/2, u, f);
        [u, p] = update(k/2, u, du);
    end

    function [u, p] = cn_step(u, f)
        du = step(k, u, f);
        [u, p] = update(k, u, du);
    end

    function [u, p] = solve(fn, u, f)
        [u, p] = feval(fn, u, f);
    end

    solver = @solve;
    handles = {@be_step, @cn_step};
end
