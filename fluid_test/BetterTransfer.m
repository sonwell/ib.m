function G = BetterTransfer(fn, N, pl, offset, F, V)
    G = feval(fn, N, pl, offset, F, V);
end

function G = l2e(N, pf, offset, vf, V)
    G = zeros(N * N * N, 1);
    O = ones(size(pf, 1), 1) * offset;
    spf = pf - O;
    indices = nearest_grid_index(N, spf);
    x = indices(:, 1);
    y = indices(:, 2);
    z = indices(:, 3);

    for i=0:63
        dx = mod(i, 4) - 1;
        dy = mod(floor(i / 4), 4) - 1;
        dz = floor(i / 16) - 1;
        
        j = 1 + mod(x + dx, N) + N * (mod(y + dy, N) + N * mod(z + dz, N));
        t = [(x + dx) / N, (y + dy) / N, (z + dz) / N];
        d = t - mod(spf, 1);
        a = accumarray(j, delta(N, d(:, 1)) .* delta(N, d(:, 2)) .* delta(N, d(:, 3)) .* vf);
        G(j) = G(j) + a(j);
    end
    G = V * G;
end

function G = e2l(N, pt, offset, vf, V)
    G = zeros(size(pt, 1), 1);
    O = ones(size(pt, 1), 1) * offset;
    spt = pt - O;
    indices = nearest_grid_index(N, spt);
    x = indices(:, 1);
    y = indices(:, 2);
    z = indices(:, 3);

    for i=0:63
        dx = mod(i, 4) - 1;
        dy = mod(floor(i / 4), 4) - 1;
        dz = floor(i / 16) - 1;

        j = 1 + mod(x + dx, N) + N * (mod(y + dy, N) + N * mod(z + dz, N));
        t = [(x + dx) / N, (y + dy) / N, (z + dz) / N];
        d = t - mod(spt, 1);
        G = G + delta(N, d(:, 1)) .* delta(N, d(:, 2)) .* delta(N, d(:, 3)) .* vf(j); 
    end
    G = V * G;
end

function d = nearest_grid_index(N, pt)
    d = floor(mod(pt, 1) * N + 1e-14);
end

function v = delta(N, d)
    v = 0.25 * N * (1 + cos(pi * d * N / 2)) .* (0.5 - abs(abs(d) - 0.5) < 2 / N);
end
