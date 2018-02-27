function G = Transfer(fn, domain, pl, offset, F)
    G = feval(fn, domain, pl, offset, F);
end

function G = l2e(domain, pf, offset, vf)
    npts = size(pf, 1);
    one = ones(npts, 1);

    nx = domain.nx;
    ny = domain.ny;
    nz = domain.nz;
    mins = one * domain.bounds(:, 1)';
    maxs = one * domain.bounds(:, 2)';
    sizes = maxs - mins;
    h = grid_size(domain);
    hs = one * h;

    G = zeros(nx * ny * nz, 1);
    O = one * (offset .* h);

    spf = mod(pf - O - mins, sizes) + mins;
    indices = nearest_grid_index(domain, spf);
    x = indices(:, 1);
    y = indices(:, 2);
    z = indices(:, 3);

    for i=0:63
        dx = mod(i, 4) - 1;
        dy = mod(floor(i / 4), 4) - 1;
        dz = floor(i / 16) - 1;
        
        j = 1 + mod(x + dx, nx) + nx * (mod(y + dy, ny) + ny * mod(z + dz, nz));
        t = [(x + dx), (y + dy), (z + dz)] .* hs;
        d = t - (spf - mins);
        pd = 0.5 * sizes - abs(abs(d) - 0.5 * sizes);
        a = accumarray(j, delta(h(1), pd(:, 1)) .* delta(h(2), pd(:, 2)) .* delta(h(3), pd(:, 3)) .* vf);
        G(j) = G(j) + a(j);
    end
end

function G = e2l(domain, pt, offset, vf)
    npts = size(pt, 1);
    one = ones(npts, 1);

    nx = domain.nx;
    ny = domain.ny;
    nz = domain.nz;
    mins = one * domain.bounds(:, 1)';
    maxs = one * domain.bounds(:, 2)';
    sizes = maxs - mins;
    h = grid_size(domain);
    hs = one * h;

    G = zeros(npts, 1);
    O = one * (offset .* h);

    spt = pt - O;
    indices = nearest_grid_index(domain, spt);
    x = indices(:, 1);
    y = indices(:, 2);
    z = indices(:, 3);

    for i=0:63
        dx = mod(i, 4) - 1;
        dy = mod(floor(i / 4), 4) - 1;
        dz = floor(i / 16) - 1;

        j = 1 + mod(x + dx, nx) + nx * (mod(y + dy, ny) + ny * mod(z + dz, nz));
        t = [(x + dx), (y + dy), (z + dz)] .* hs;
        d = t - mod(spt - mins, sizes);
        pd = 0.5 * sizes - abs(abs(d) - 0.5 * sizes);
        G = G + delta(h(1), pd(:, 1)) .* delta(h(2), pd(:, 2)) .* delta(h(3), pd(:, 3)) .* vf(j); 
    end
end

function h = grid_size(domain)
    nx = domain.nx;
    ny = domain.ny;
    nz = domain.nz;

    mins = domain.bounds(:, 1)';
    maxs = domain.bounds(:, 2)';
    sizes = maxs - mins;
    h = sizes ./ [nx ny nz];
end

function d = nearest_grid_index(domain, pt)
    npts = size(pt, 1);
    one = ones(npts, 1);
    h = one * grid_size(domain);
    mins = one * domain.bounds(:, 1)';
    d = floor((pt - mins) ./ h);
end

function v = delta(h, d)
    v = 0.25 / h * (1 + cos(pi * d / (2 * h))) .* (d < 2 * h);
end
