classdef PeriodicCylinder < ClosedSurface
    methods
        function obj = PeriodicCylinder(n, m)
            rbf = PolyharmonicSpline(4);
            obj@ClosedSurface(n, m, rbf);
            obj.ds = 2 * pi * obj.ds;
        end

        function [id, da, db, daa, dab, dbb, phi, psi] = operators(obj, rbf, data, sample)
            npts = size(data, 1);
            r = obj.metric(data, data);
            phi = rbf.phi(r);
            one = ones(npts, 1);
            itp = [phi one data(:, 2); one' 0 0; data(:, 2)' 0 0];
            trim = @(m) m(:, 1:npts);

            mpts = size(sample, 1);
            [r, r_a, r_b, r_aa, r_ab, r_bb] = obj.metric(data, sample);
            psi = rbf.phi(r);
            one = ones(mpts, 1);
            zero = zeros(mpts, 1);
            id = trim((itp' \ [psi one sample(:, 2)]')');
            da = trim((itp' \ [rbf.dphi(r, r_a) zero zero]')');
            db = trim((itp' \ [rbf.dphi(r, r_b) zero one]')');
            daa = trim((itp' \ [rbf.ddphi(r, r_aa, r_a .* r_a) zero zero]')');
            dab = trim((itp' \ [rbf.ddphi(r, r_ab, r_a .* r_b) zero zero]')');
            dbb = trim((itp' \ [rbf.ddphi(r, r_bb, r_b .* r_b) zero zero]')');
        end
    end

    methods(Static)
        function [r, r_a, r_b, r_aa, r_ab, r_bb] = metric(data, sample)
            [dt, dz] = process(data, sample);

            r = sqrt(2*(1-cos(dt)) + 2 * (1-cos(2*pi*dz)));
            r_a = sin(dt);
            r_b = 2 * pi * sin(2*pi*dz);
            r_aa = cos(dt);
            r_ab = 0;
            r_bb = 4*pi^2*cos(2*pi*dz);
        end

        function x = shape(params)
            theta = params(:, 1);
            z = params(:, 2);

            x = [cos(theta) sin(theta) z];
        end

        function params = sample(n)
            sc = 2 * pi;
            layers = ceil(sqrt(n / sc));
            uniform = ((1:n)-1)' / n;
            z = uniform;
            theta = mod(2 * pi * layers * uniform, 2 * pi);
            params = [theta z];
        end
    end
end

function [dt, dz] = process(data, sample)
    alpha = data(:, 1);
    zeta = data(:, 2);
    gamma = sample(:, 1);
    eta = sample(:, 2);

    [tj, tk] = ndgrid(gamma, alpha);
    [zj, zk] = ndgrid(eta, zeta);
    dt = tj - tk;
    dz = zj - zk;
end
