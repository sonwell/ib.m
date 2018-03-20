classdef RedBloodCell < Sphere & ForceMixin
    methods
        function obj = RedBloodCell(n, m, varargin)
            data = RedBloodCell.sample(n);
            sample = RedBloodCell.sample(m);
            obj@Sphere(data, sample);
            obj@ForceMixin(data, varargin{:});

            x = RedBloodCell.shape(data);
            y = Sphere.shape(data);
            gx = obj.geometry(x);
            gy = obj.geometry(y);
            obj.ds = sqrt((gx.E .* gx.G - gx.F .^ 2) ./ (gy.E .* gy.G - gy.F .^ 2)).* obj.ds;
        end
    end

    methods(Static)
        function x = shape(params)
            xs = Sphere.shape(params);

            R = 4e-4;  % radius of 4e-4 cm
            r = xs(:, 1).^2 + xs(:, 3).^2;
            xt = R * xs(:, 1);
            % Coefficients from Skalak, et al., (1973).
            yt = R / 2 * xs(:, 2) .* (0.21 + 2.0 * r - 1.12 * r.^2);
            zt = R * xs(:, 3);
            x = [2^-9 + xt, 193/350 * R +  yt, 2^-9 + zt];
        end
    end
end
