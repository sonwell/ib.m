classdef Torus < ClosedSurface
    methods (Static)
        function params = sample(n)
            layers = ceil(sqrt(n));
            uniform = ((1:n)-1)' / n;
            alpha = mod(2 * pi * layers * uniform, 2 * pi);
            beta = mod(2 * pi * uniform, 2 * pi);
            params = [alpha beta];
        end

        function [r, rdr, drdr] = metric(data, sample)
            ad = data(:, 1);
            bd = data(:, 2);
            as = sample(:, 1);
            bs = sample(:, 2);

            [asr, adr] = ndgrid(as, ad);
            [bsr, bdr] = ndgrid(bs, bd);
            da = asr - adr;
            db = bsr - bdr;

            r = sqrt(2*(1-cos(da)) + 2 * (1-cos(db)));
            r_a = sin(da);
            r_b = sin(db);
            r_aa = cos(da);
            r_ab = 0;
            r_bb = cos(db);

            rdr = {r_a, r_b};
            drdr = {r_aa, r_ab; r_ab, r_bb};
        end
    end

    methods
        function self = Torus(n, m, rbf, poly, varargin)
            qpoly = Polynomials(0, poly.dimensions);
            self@ClosedSurface(n, m, rbf, poly, qpoly, 4*pi^2, @(d) 1, varargin{:});
        end

        function x = shape(~, params)
            theta = params(:, 1);
            phi = params(:, 2);
            x = [(1 + cos(theta)).*cos(phi) (1 + cos(theta)).*sin(phi) sin(theta)];
        end
    end
end
