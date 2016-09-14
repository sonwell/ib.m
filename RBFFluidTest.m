function RBFFluidTest(N, CP, EP, rho, mu)
    alpha = CP(:, 1);
    beta = CP(:, 2);
    gamma = EP(:, 1);
    delta = EP(:, 2);

    R = 0.25;
    x = sin(beta) .* cos(alpha);
    z = cos(beta);
    r = x .^ 2 + z .^ 2;
    y = 0.5 * sin(beta) .* sin(alpha) .* (0.32 + 2.003 * r - 1.123 * r .^ 2);

    x = 0.5 + R * x;
    y = 0.5 + R * y;
    z = 0.5 + R * z;

    y0 = y;

    o = ones(N, 1);
    i = (0:(N-1))' / N;
    O = kron(o, kron(o, o));
    c = 0.5 / N + [kron(o, kron(o, i)), kron(o, kron(i, o)), kron(i, kron(o, o))];
    gx = c + O * [-0.5 / N 0 0];
    gy = c + O * [0 -0.5 / N 0];
    gz = c + O * [0 0 -0.5 / N];

    Force = ForceCalculator(alpha, beta, gamma, delta, x, y, z, 1e-1, 1e-2, 3.9e-7);
    pts = size(EP, 1);
    sa = 4 * pi * R^2 / pts;

    Ue = zeros(N^3, 3);
    Ff = zeros(N^3, 3);

    sps = 1000000;
    k = 1 / sps;
    eps = 0.001;
    Z = gy(:, 3);
    cell = fopen('cell3.txt', 'w');
    fluid = fopen('fluid3.txt', 'w');

    Ue(:, 2) = 1; % 0.1 * sin(2 * pi * Z);
    Ff(:, 2) = 0; % 0.1 * 4 * pi^2 * sin(2 * pi * Z);
    Solver = StokesSolver(rho, mu, k, N);
    for m = 1:30
        for n = 1:sps
            [Xl, Fls, Flb] = Force(x, y, z);
            Fl = Fls + Flb;

            Fex = BetterTransfer('l2e', N, Xl, [0 0.5 0.5] / N, Fl(:, 1), sa);
            Fey = BetterTransfer('l2e', N, Xl, [0.5 0 0.5] / N, Fl(:, 2), sa);
            Fez = BetterTransfer('l2e', N, Xl, [0.5 0.5 0] / N, Fl(:, 3), sa);
            Fe = [Fex Fey Fez];

            [Us, Ps] = Solver('be', Ue, Fe + Ff);
            Ulx = BetterTransfer('e2l', N, [x y z], [0 0.5 0.5] / N, Us(:, 1), 1 / N^3);
            Uly = BetterTransfer('e2l', N, [x y z], [0.5 0 0.5] / N, Us(:, 2), 1 / N^3);
            Ulz = BetterTransfer('e2l', N, [x y z], [0.5 0.5 0] / N, Us(:, 3), 1 / N^3);

            xi = x + k / 2 * Ulx;
            yi = y + k / 2 * Uly;
            zi = z + k / 2 * Ulz;

            [Xl, Fls, Flb, curr, H] = Force(xi, yi, zi);
            Fl = Fls + Flb;

            Fex = BetterTransfer('l2e', N, Xl, [0 0.5 0.5] / N, Fl(:, 1), sa);
            Fey = BetterTransfer('l2e', N, Xl, [0.5 0 0.5] / N, Fl(:, 2), sa);
            Fez = BetterTransfer('l2e', N, Xl, [0.5 0.5 0] / N, Fl(:, 3), sa);
            Fe = [Fex Fey Fez];

            [Ue, Pe] = Solver('cn', Ue, Fe + Ff);
            Ulx = BetterTransfer('e2l', N, [x y z], [0 0.5 0.5] / N, Ue(:, 1), 1 / N^3);
            Uly = BetterTransfer('e2l', N, [x y z], [0.5 0 0.5] / N, Ue(:, 2), 1 / N^3);
            Ulz = BetterTransfer('e2l', N, [x y z], [0.5 0.5 0] / N, Ue(:, 3), 1 / N^3);

            x = x + k * Ulx;
            y = y + k * Uly;
            z = z + k * Ulz;

%            data = [((m - 1) * sps + n) * ones(pts, 1) (1:pts)' x y z Ulx Uly Ulz Fls Flb];
%            fprintf(cell, '%d %d %0.15f %0.15f %0.15f %0.15f %0.15f %0.15f %0.15f %0.15f %0.15f %0.15f %0.15f %0.15f\n', data');
%            data = [((m - 1) * sps + n) * ones(N * N * N, 1) (1:(N * N * N))' c Ue Pe];
%            fprintf(fluid, '%d %d %0.15f %0.15f %0.15f %0.15f %0.15f %0.15f %0.15f\n', data');
            scatter3(Xl(:, 1), Xl(:, 2), Xl(:, 3), 200, H, 'filled');
            axis equal; axis vis3d;
            drawnow;
        end
    end
end
