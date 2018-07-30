classdef PeriodicSheet < Torus
    methods
        function obj = PeriodicSheet(n, m, varargin)
            rbf = PolyharmonicSpline(4);
            poly = Polynomials(1, 2);
            obj@Torus(n, m, rbf, poly, varargin{:});
        end

        function x = shape(~, params)
            y = params(:, 1) / (2 * pi);
            z = params(:, 2) / (2 * pi);

            x = [y z 0*y];
        end
    end
end
