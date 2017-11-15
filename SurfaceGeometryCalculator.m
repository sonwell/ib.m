function calculator = SurfaceGeometryCalculator(Da, Db, Daa, Dab, Dbb)
    function S = geometry(x, y, z)
        X = [x y z];
        ta = Da * X;
        tb = Db * X;

        taa = Daa * X;
        tab = Dab * X;
        tbb = Dbb * X;

        m = cross(ta, tb, 2);

        E = dot(ta, ta, 2);
        F = dot(ta, tb, 2);
        G = dot(tb, tb, 2);

        J = sqrt(E .* G - F .^ 2);
        n = [1./J 1./J 1./J] .* m;

        e = dot(taa, n, 2);
        f = dot(tab, n, 2);
        g = dot(tbb, n, 2);

        H = (e .* G - 2 * f .* F + g .* E) ./ (2 * (E .* G - F .^ 2));
        K = (e .* g - f .^ 2) ./ (E .* G - F .^ 2);

        S.ta = ta;
        S.tb = tb;
        S.n = n;

        S.taa = taa;
        S.tab = tab;
        S.tbb = tbb;

        S.E = E; S.e = e;
        S.F = F; S.f = f;
        S.G = G; S.g = g;
        S.H = H;
        S.K = K;

        S.Ea = 2 * dot(taa, ta, 2);
        S.Eb = 2 * dot(tab, ta, 2);
        S.Fa = dot(taa, tb, 2) + dot(ta, tab, 2);
        S.Fb = dot(tab, tb, 2) + dot(ta, tbb, 2);
        S.Ga = 2 * dot(tab, tb, 2);
        S.Gb = 2 * dot(tbb, tb, 2);
    end
    calculator = @geometry;
end
