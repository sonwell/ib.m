classdef Sphere < ClosedSurface
    methods
        function obj = Sphere(data, sample)
            rbf = PolyharmonicSpline(20);

            npts = size(data, 1);
            r = Sphere.metric(data, data);
            phi = rbf.phi(r);
            one = ones(npts, 1);
            zero = zeros(npts, 1);
            itp = [phi one; one' 0];
            trim = @(m) m(:, 1:npts);
            ds1 = itp' \ [zero; 4*pi];

            mpts = size(sample, 1);
            [r, r_a, r_b, r_aa, r_ab, r_bb] = Sphere.metric(data, sample);
            one = ones(mpts, 1);
            zero = zeros(mpts, 1);
            id = trim((itp' \ [rbf.phi(r) one]')');
            da = trim((itp' \ [rbf.dphi(r, r_a) zero]')');
            db = trim((itp' \ [rbf.dphi(r, r_b) zero]')');
            daa = trim((itp' \ [rbf.ddphi(r, r_aa, r_a .* r_a) zero]')');
            dab = trim((itp' \ [rbf.ddphi(r, r_ab, r_a .* r_b) zero]')');
            dbb = trim((itp' \ [rbf.ddphi(r, r_bb, r_b .* r_b) zero]')');
            ds = id * (spdiags(1 ./ sum(id)', [0], npts, npts) * ds1(1:npts));

            obj@ClosedSurface(data, sample, id, da, db, daa, dab, dbb, phi, ds(1:mpts));
        end
    end

    methods(Static)
        function [r, r_a, r_b, r_aa, r_ab, r_bb] = metric(data, sample)
            [ad, bj, bk] = process(data, sample);

            r    = sqrt(2 * (1-sin(bj).*sin(bk).*cos(ad)-cos(bj).*cos(bk)) + 1e-15);
            r_a  = sin(bj).*sin(bk).*sin(ad);
            r_b  = sin(bj).*cos(bk) -cos(bj).*sin(bk).*cos(ad);
            r_aa = sin(bj).*sin(bk).*cos(ad);
            r_ab = cos(bj).*sin(bk).*sin(ad);
            r_bb = cos(bj).*cos(bk) +sin(bj).*sin(bk).*cos(ad);
        end

        function x = shape(params)
            theta = params(:, 1);
            phi = params(:, 2);

            x = [sin(phi) .* cos(theta), sin(phi) .* sin(theta), cos(phi)];
        end

        function params = sample(n)
            rotations = sqrt(n);
            uniform = ((1:n)' - 0.5) / n;
            spacing = sqrt(csc(pi * uniform));
            scaling = 1 / sum(spacing);
            phi = pi * (cumsum(scaling * spacing) - 0.5 * scaling * spacing);
            theta = mod(2 * rotations * phi, 2 * pi);
            params = [theta phi];
        end
    end 
end

function [da, bj, bk] = process(data, sample)
    alpha = data(:, 1);
    beta = data(:, 2);
    gamma = sample(:, 1);
    delta = sample(:, 2);

    [aj, ak] = ndgrid(gamma, alpha);
    [bj, bk] = ndgrid(delta, beta);
    da = aj - ak;
end
