classdef ClosedSurface < ParamObject
    % ClosedSurface Geometric representation of a closed shape in 3D.

    properties
        data_sites  % Interpolation point parameters
        sample_sites  % Evaluation point parameters
        rbf
        poly
        id  % Evaluation matrix
        da  % Differentiation matrix with respect to the first parameter
        db  % Differentiation matrix with respect to the second parameter
        daa  % Second derivative matrix with respect to the first parameter
        dab  % Mixed second derivative matrix
        dbb  % Second derivative matrix with respect to the second parameter
        ds  % "Element" surface area for each evaluation point
        phi
        psi
    end

    methods (Abstract)
        shape(obj, x)
    end

    methods (Abstract, Static)
        sample(n)
        metric(data, sample)
    end

    methods
        function obj = ClosedSurface(n, m, rbf, poly, varargin)
            obj@ParamObject(varargin{:});
            data = obj.sample(n);
            sample = obj.sample(m);
            [id, da, db, daa, dab, dbb, phi, psi] = obj.operators(rbf, poly, data, sample);

            obj.rbf = rbf;
            obj.poly = poly;
            obj.data_sites = data;
            obj.sample_sites = sample;
            obj.id = id;
            obj.da = da;
            obj.db = db;
            obj.daa = daa;
            obj.dab = dab;
            obj.dbb = dbb;
            obj.phi = phi;
            obj.psi = psi;

            obj.ds = obj.surface_area(phi, psi);
        end

        function sa = surface_area(obj, phi, psi)
            n = size(phi, 1);
            m = size(psi, 1);
            trim = @(M) M(:, 1:n);
            itp = [phi ones(n, 1); ones(1, n) 0];
            id = trim((itp' \ [psi ones(m, 1)]')');
            sc = spdiags(1 ./ sum(id)', [0], n, n);
            w = itp' \ [zeros(n, 1); 1];
            sa = id * (sc * w(1:n));
        end

        function [px, py, pz] = surf(obj, x)
            data = obj.data_sites;
            rbf = obj.rbf;
            n = size(data, 1);
            phi = obj.phi;
            one = ones(n, 1);
            itp = [phi one; one' 0];
            w = itp \ [x; zeros(1, size(x, 2))];

            function y = surfer(u, v, c)
                m = numel(u);
                r = obj.metric(data, [u(:) v(:)]);
                psi = rbf.phi(r);
                one = ones(m, 1);
                y = [psi one] * c;
            end
            
            px = @(u, v) reshape(surfer(u, v, w(:, 1)), size(u));
            py = @(u, v) reshape(surfer(u, v, w(:, 2)), size(u));
            pz = @(u, v) reshape(surfer(u, v, w(:, 3)), size(u));
        end

        function r = measure_ratio(obj, x, y)
            gx = obj.geometry(x);
            gy = obj.geometry(y);
            r = sqrt(gx.I ./ gy.I);
        end

        function [x, ta, tb, taa, tab, tbb] = representation(obj, y)
            x = obj.id * y;

            ta = obj.da * y;
            tb = obj.db * y;

            taa = obj.daa * y;
            tab = obj.dab * y;
            tbb = obj.dbb * y;
        end

        function geom = geometry(obj, y)
            [x, ta, tb, taa, tab, tbb] = obj.representation(y);
            geom = obj.parameters(ta, tb, taa, tab, tbb);
            geom.x = x;
        end

        function [id, da, db, daa, dab, dbb, phi, psi] = operators(obj, rbf, poly, data, sample)
            npts = size(data, 1);
            r = obj.metric(data, data);
            phi = rbf.phi(r);
            p = poly.p(data);
            pad = zeros(size(p, 2));
            itp = [phi p; p' pad];
            trim = @(m) m(:, 1:npts);

            [r, r_a, r_b, r_aa, r_ab, r_bb] = obj.metric(data, sample);
            psi = rbf.phi(r);
            q = poly.p(sample);
            qa = poly.dp(sample, 1);
            qb = poly.dp(sample, 2);
            qaa = poly.dp(sample, [1 1]);
            qab = poly.dp(sample, [1 2]);
            qbb = poly.dp(sample, [2 2]);

            id = trim((itp' \ [psi q]')');
            da = trim((itp' \ [rbf.dphi(r, r_a) qa]')');
            db = trim((itp' \ [rbf.dphi(r, r_b) qb]')');
            daa = trim((itp' \ [rbf.ddphi(r, r_aa, r_a .* r_a) qaa]')');
            dab = trim((itp' \ [rbf.ddphi(r, r_ab, r_a .* r_b) qab]')');
            dbb = trim((itp' \ [rbf.ddphi(r, r_bb, r_b .* r_b) qbb]')');
        end
    end

    methods (Static)
        function geom = parameters(ta, tb, taa, tab, tbb)
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

            geom.I = I;
            geom.II = II;
            geom.H = H;
            geom.K = K;

            geom.Ea = 2 * dot(taa, ta, 2);
            geom.Eb = 2 * dot(tab, ta, 2);
            geom.Fa = dot(taa, tb, 2) + dot(ta, tab, 2);
            geom.Fb = dot(tab, tb, 2) + dot(ta, tbb, 2);
            geom.Ga = 2 * dot(tab, tb, 2);
            geom.Gb = 2 * dot(tbb, tb, 2);
        end
    end
end
