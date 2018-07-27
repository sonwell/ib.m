classdef PeriodicCylinder < Torus
    methods
        function obj = PeriodicCylinder(n, m, varargin)
            rbf = PolyharmonicSpline(4);
            poly = PolynomialSubset(1, 2, 2);
            obj@Torus(n, m, rbf, poly, varargin{:});
        end

        function x = shape(~, params)
            theta = params(:, 1);
            z = params(:, 2) / (2 * pi);
            x = [cos(theta) sin(theta) z];
        end
    end
end
