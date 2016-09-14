function solver = StokesSolver(rho, mu, k, N)
    G = spdiags(ones(N, 2), [-N + 1, 1], N, N);
    E = spdiags(ones(N, 4), [-N+1, -1, 1, N-1], N, N);
    A = spdiags(ones(N, 1) * [1 1 -6 1 1], [-N+1, -1, 0, 1, N-1], N, N);
    I = speye(N);
    I3 = kron(I, kron(I, I));
    L = N * N * (kron(E, kron(I, I)) + kron(I, kron(E, I)) + kron(I, kron(I, A)));
    
    lambda = mu * k / (2 * rho);
    Hp = I3 + lambda * L;
    Hm = I3 - lambda * L;

    %Lf = chol(-L);
    %Hf = chol(Hm);

    Dx = N * (kron(I, kron(I, G)) - I3);
    Dy = N * (kron(I, kron(G, I)) - I3);
    Dz = N * (kron(G, kron(I, I)) - I3);
    
    function [U, P] = solve(fn, U, F)
        U = feval(fn, rho, Hp, Hm, k, U, F);
        dU = Dx * U(:, 1) + Dy * U(:, 2) + Dz * U(:, 3);
        dU = dU - mean(dU);
        dtphi = L \ dU;
        U = U + [Dx' * dtphi, Dy' * dtphi, Dz' * dtphi];
        P = Hm * dtphi / k;
    end

    solver = @solve;
end

function U = be(rho, Hp, Hm, k, U, F)
    f = k / (2 * rho) * F;
    U = Hm \ (U + f);
end

function U = cn(rho, Hp, Hm, k, U, F)
    f = k / rho * F;
    U = Hm \ (Hp * U + f);
end
