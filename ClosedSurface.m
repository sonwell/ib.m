classdef ClosedSurface
    properties
        data_sites
        sample_sites
        id
        da
        db
        daa
        dab
        dbb
        phi
        ds
    end

    methods
        function obj = ClosedSurface(data, sample, id, da, db, daa, dab, dbb, phi, ds)
            obj.data_sites = data;
            obj.sample_sites = sample;
            obj.id = id;
            obj.da = da;
            obj.db = db;
            obj.daa = daa;
            obj.dab = dab;
            obj.dbb = dbb;
            obj.phi = phi;
            obj.ds = ds;
        end

        function geom = parameters(obj, ta, tb, taa, tab, tbb)
            E = dot(ta, ta, 2);
            F = dot(ta, tb, 2);
            G = dot(tb, tb, 2);

            I = E .* G - F .^ 2;
            J = sqrt(I);
            n = cross(ta, tb, 2) ./ [J J J];

            e = dot(taa, n, 2);
            f = dot(tab, n, 2);
            g = dot(tbb, n, 2);

            II = e .* g - f .^ 2;
            H = (e .* G - 2 * f .* F + g .* E) ./ (2 * I);
            K = II ./ I;

            geom.ta = ta;
            geom.tb = tb;
            geom.n = n;

            geom.taa = taa;
            geom.tab = tab;
            geom.tbb = tbb;

            geom.E = E;
            geom.F = F;
            geom.G = G;

            geom.e = e;
            geom.f = f;
            geom.g = g;

            geom.H = H;
            geom.K = K;

            geom.Ea = 2 * dot(taa, ta, 2);
            geom.Eb = 2 * dot(tab, ta, 2);
            geom.Fa = dot(taa, tb, 2) + dot(ta, tab, 2);
            geom.Fb = dot(tab, tb, 2) + dot(ta, tbb, 2);
            geom.Ga = 2 * dot(tab, tb, 2);
            geom.Gb = 2 * dot(tbb, tb, 2);
        end
    
        function [x, ta, tb, taa, tab, tbb] = representation(obj, y)
            x = obj.id * y;

            ta = obj.da * y;
            tb = obj.db * y;

            taa = obj.daa * y;
            tab = obj.dab * y;
            tbb = obj.dbb * y;
        end

        function geom = geometry(obj, x)
            [x, ta, tb, taa, tab, tbb] = obj.representation(x);
            geom = obj.parameters(ta, tb, taa, tab, tbb);
            geom.x = x;
        end
    end
end
