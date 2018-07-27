classdef RedBloodCell < Sphere & ForceMixin
    methods
        function obj = RedBloodCell(n, m, varargin)
            obj@Sphere(n, m);
            obj@ForceMixin(varargin{:});
        end

        function x = shape(obj, params)
            xs = shape@Sphere(obj, params);

            R = 3.91e-4;  % radius of 3.91e-4 cm
            r = xs(:, 1).^2 + xs(:, 3).^2;
            xt = R * xs(:, 1);
            % Coefficients from Skalak, et al., (1973).
            yt = R / 2 * xs(:, 2) .* (0.21 + 2.0 * r - 1.12 * r.^2);
            zt = R * xs(:, 3);
            x = [2^-9 + xt, 193/350 * R +  yt, 2^-9 + zt];
        end
    end
end
