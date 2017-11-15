function [U, P] = Stokes(fn, rho, mu, k, N, U, F)
    G = spdiags(ones(N, 2), [-N + 1, 1], N, N);
    E = spdiags(ones(N, 4), [-N+1, -1, 1, N-1], N, N);
    A = spdiags(ones(N, 1) * [1 1 -6 1 1], [-N+1, -1, 0, 1, N-1], N, N);
    I = speye(N);
    I3 = kron(I, kron(I, I));
    L = N * N * (kron(E, kron(I, I)) + kron(I, kron(E, I)) + kron(I, kron(I, A)));

    [U, Hm] = feval(fn, rho, mu, k, I3, L, U, F);

    Dx = N * (kron(I, kron(I, G)) - I3);
    Dy = N * (kron(I, kron(G, I)) - I3);
    Dz = N * (kron(G, kron(I, I)) - I3);

    dU = Dx * U(:, 1) + Dy * U(:, 2) + Dz * U(:, 3);
    dU = dU - mean(dU);
    dtphi = L \ dU;
    U = U + [Dx' * dtphi, Dy' * dtphi, Dz' * dtphi];
    P = Hm * dtphi / k;
end

function [U, Hm] = be(rho, mu, k, I, L, U, F)
    lambda = mu / rho * k / 2;
    Hp = I;
    Hm = I - lambda * L;

    f = k / rho * F / 2;
    U = Hm \ (Hp * U + f);
end

function [U, Hm] = cn(rho, mu, k, I, L, U, F)
    lambda = mu / rho * k / 2;
    Hp = I + lambda * L;
    Hm = I - lambda * L;

    f = k / rho * F;
    U = Hm \ (Hp * U + f);
end
