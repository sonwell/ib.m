function X = RBCShape(params)
    alpha = params(:, 1);
    beta = params(:, 2);
    R = 2^-11;
    xs = sin(beta) .* cos(alpha);
    ys = sin(beta) .* sin(alpha);
    zs = cos(beta);

    r = (xs.^2 + zs.^2);
    xt = 2^-9 + R * xs;
    yt = 2^-9 + R / 2 * ys .* (0.21 + 2.0 * r - 1.12 * r.^2);
    zt = 2^-9 + R * zs;
    X = [xt, yt, zt];
end
