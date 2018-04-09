function [solver, L, handles] = NavierStokesSolver(domain, rho, mu, k)
    n = [domain.nx domain.ny domain.nz];
    h = (domain.bounds(:, 2) - domain.bounds(:, 1))' ./ n;

    multiplex = @(fn, n, h) deal(fn(n(1), h(1)), fn(n(2), h(2)), fn(n(3), h(3)));
    identity = @(n, h) speye(n);
    differential = @(n, h) spdiags(ones(n, 1) * [1 -1 1] / h, [-n+1, 0, 1], n, n);
    laplacian = @(n, h) spdiags(ones(n, 1) * [1 1 -2 1 1] / h^2, [-n+1, -1, 0, 1, n-1], n, n);
    average = @(n, h) spdiags(ones(n, 1) * [0.5 0.5 0.5], [-1, 1, n-1], n, n);

    [Gx, Gy, Gz] = multiplex(differential, n, h);
    [Lx, Ly, Lz] = multiplex(laplacian, n, h);
    [Ix, Iy, Iz] = multiplex(identity, n, h);
    [Ax, Ay, Az] = multiplex(average, n, h);

    I = kron(Iz, kron(Iy, Ix));
    L = kron(Lz, kron(Iy, Ix)) + kron(Iz, kron(Ly, Ix)) + kron(Iz, kron(Iy, Lx));
    Dx = kron(Iz, kron(Iy, Gx));
    Dy = kron(Iz, kron(Gy, Ix));
    Dz = kron(Gz, kron(Iy, Ix));
    Mx = kron(Iz, kron(Iy, Ax));
    My = kron(Iz, kron(Ay, Ix));
    Mz = kron(Az, kron(Iy, Ix));
    nL = -L;
    
    lambda = mu * k / (2 * rho);
    H = I - lambda * L;
    C = ichol(H);
    C2 = ichol(nL);

    function udu = momentum(ul)
        u = ul(:, 1);
        v = ul(:, 2);
        w = ul(:, 3);

        axu = Ax * u; axv = Ax * v; axw = Ax * w;
        ayu = Ay * u; ayv = Ay * v; ayw = Ay * w;
        azu = Az * u; azv = Az * v; azw = Az * w;

        uu = axu .* axu;
        uv = ayu .* axv;
        uw = azu .* axw;
        vv = ayv .* ayv;
        vw = azv .* ayw;
        ww = azw .* azw;

        hu = Dx * (uu) + Dy * (uv) + Dz * (uw);
        hv = Dx * (uv) + Dy * (vv) + Dz * (vw);
        hw = Dx * (uw) + Dy * (vw) + Dz * (ww);

        udu = [hu hv hw];
    end

    function du = step(k, u, ul, f)
        h = momentum(ul);
        rhs = k / rho * (mu * (L * u) + h + f);

        du = zeros(size(u));
        [dux, flx, rr, itx, rv] = pcg(H, rhs(:, 1), 1e-10, 100000, C, C');
        [duy, fly, rr, ity, rv] = pcg(H, rhs(:, 2), 1e-10, 100000, C, C');
        [duz, flz, rr, itz, rv] = pcg(H, rhs(:, 3), 1e-10, 100000, C, C');
        du(:, 1) = dux;
        du(:, 2) = duy;
        du(:, 3) = duz;
    end

    function [u, p] = update(k, u, du)
        u_star = u + du;
        div_u_star = Dx * u_star(:, 1) + Dy * u_star(:, 2) + Dz * u_star(:, 3);
        pr = div_u_star - mean(div_u_star);
        [dtphi, fl, rr, it, rv] = pcg(nL, -pr, 1e-10, 100000, C2, C2');
        u = u_star + [Dx' * dtphi, Dy' * dtphi, Dz' * dtphi];
        p = H * dtphi / k;
    end

    function [u, p] = be_step(u, ul, f)
        du = step(k/2, u, f);
        [u, p] = update(k/2, u, du);
    end

    function [u, p] = cn_step(u, ul, f)
        du = step(k, u, ul, f);
        [u, p] = update(k, u, du);
    end

    function [u, p] = solve(fn, u, ul, f)
        [u, p] = feval(fn, u, ul, f);
    end

    solver = @solve;
    handles = {@be_step, @cn_step};
end
