function RBFFluidTest(N, CP, EP, rho, mu)
    alpha = CP(:, 1);
    beta = CP(:, 2);
    gamma = EP(:, 1);
    delta = EP(:, 2);

    % initial configuration
    R = 0.25;
    x = sin(beta) .* cos(alpha);
    y = sin(beta) .* sin(alpha);
    z = cos(beta);

    x = 0.5 + R * x;
    y = 0.5 + R * y;
    z = 0.5 + R * z;

    o = ones(N, 1);
    i = (0:(N-1))' / N;
    O = kron(o, kron(o, o));
    c = 0.5 / N + [kron(o, kron(o, i)), kron(o, kron(i, o)), kron(i, kron(o, o))];
    gx = c + O * [-0.5 / N 0 0];
    gy = c + O * [0 -0.5 / N 0];
    gz = c + O * [0 0 -0.5 / N];

    Force = ForceCalculator(alpha, beta, gamma, delta, x, y, z, 1e-3, 1e-2, 3.9e-7);
    pts = size(EP, 1);
    sa = 4 * pi * R^2 / pts;

    Ue = zeros(N^3, 3);
    Ff = zeros(N^3, 3);

    sps = 1000000;
    k = 1 / sps;
    eps = 0.001;
    Z = gy(:, 3);

    Ue(:, 2) = 1;
    Ff(:, 2) = 0;
    Solver = StokesSolver(rho, mu, k, N);
    for m = 1:30
        for n = 1:sps
            [Xl, Fls, Flb, curr] = Force(x, y, z);
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

            [Xl, Fls, Flb] = Force(xi, yi, zi);
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

            scatter3(Xl(:, 1), Xl(:, 2), Xl(:, 3), 200, sin(gamma) .* sin(delta), 'filled')
            axis equal; axis vis3d;
            colorbar;
            drawnow;
        end
    end
end
