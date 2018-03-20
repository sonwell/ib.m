classdef BloodVessel < PeriodicCylinder & ForceMixin
    methods
        function obj = BloodVessel(n, m, varargin)
            data = BloodVessel.sample(n);
            sample = BloodVessel.sample(m);
            obj@PeriodicCylinder(data, sample);
            obj@ForceMixin(data, varargin{:});

            x = BloodVessel.shape(data);
            y = PeriodicCylinder.shape(data);
            gx = obj.geometry(x);
            gy = obj.geometry(y);
            obj.ds = sqrt((gx.E .* gx.G - gx.F.^2) ./ (gy.E .* gy.G - gy.F.^2)) .* obj.ds;
        end    
    end

    methods(Static)
        function x = shape(params)
            theta = params(:, 1);
            zeta = params(:, 2);

            rho = @(x) tanh(10*(x+0.5))-tanh(10*(x-0.5));
            r = 0;
            for n = -1:1
                r = r + rho(-1 + 2 * (zeta - n));
            end

            x = [2^-9 + 9e-4 * (1 - 0.25 * r) .* cos(theta), ...
                 2^-8 * zeta, ...
                 2^-9 + 9e-4 * (1 - 0.25 * r) .* sin(theta)];
        end

        function params = sample(n)
            rho = @(x) tanh(10*(x+0.5))-tanh(10*(x-0.5));
            drho = @(x) -tanh(10*(x+0.5)).^2+tanh(10*(x-0.5)).^2;
            r = 0;
            dr = 0;

            sc = 2 * pi / 7.8125;
            layers = floor(sqrt(n / sc));
            uniform = ((1:n)-0.5)' / n;
            for n = -1:1
                r = r + rho(-1 + 2 * (uniform - n));
                dr = dr + drho(-1 + 2 * (uniform - n));
            end
            spacing = ((1-0.25 * r).^2 .* ((1-0.25 * dr).^2+1)).^(-1/4);
            scaling = 1 / sum(spacing);
            z = cumsum(scaling * spacing) - scaling * spacing;
            theta = mod(2 * pi * layers * z, 2 * pi);
            params = [theta z];
        end
    end
end
