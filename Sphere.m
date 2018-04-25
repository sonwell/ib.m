classdef Sphere < ClosedSurface
    methods
        function obj = Sphere(n, m, varargin)
            rbf = PolyharmonicSpline(10);
            poly = Polynomials(0, 2);
            obj@ClosedSurface(n, m, rbf, poly, varargin{:});
            obj.ds = 4 * pi * obj.ds;
        end

        function [x, r] = shape(~, params)
            phi = params(:, 1);
            theta = params(:, 2);

            x = [sin(phi) .* cos(theta), sin(phi) .* sin(theta), cos(phi)];
            r = 1;
        end
    end

    methods(Static)
        function [r, r_p, r_t, r_pp, r_pt, r_tt] = metric(data, sample)
            [pj, pk, dt] = process(data, sample);

            r    = sqrt(2 * (1-sin(pj).*sin(pk).*cos(dt)-cos(pj).*cos(pk)) + 2^-51);
            r_p  = sin(pj).*cos(pk) -cos(pj).*sin(pk).*cos(dt);
            r_t  = sin(pj).*sin(pk).*sin(dt);
            r_pp = cos(pj).*cos(pk) +sin(pj).*sin(pk).*cos(dt);
            r_pt = cos(pj).*sin(pk).*sin(dt);
            r_tt = sin(pj).*sin(pk).*cos(dt);
        end

        function params = sample(n)
            rotations = sqrt(n);
            uniform = linspace(1/(2*n), 1-1/(2*n), n)';
            spacing = sqrt(csc(pi * uniform));
            scaling = 1 / sum(spacing);
            phi = pi * (cumsum(scaling * spacing) - 0.5 * scaling * spacing);
            theta = mod(2 * rotations * phi, 2 * pi);
            params = [phi theta];
        end
    end 
end

function [pj, pk, dt] = process(data, sample)
    phi_d = data(:, 1);
    theta_d = data(:, 2);
    phi_s = sample(:, 1);
    theta_s = sample(:, 2);

    [pj, pk] = ndgrid(phi_s, phi_d);
    [tj, tk] = ndgrid(theta_s, theta_d);
    dt = tj - tk;
end
