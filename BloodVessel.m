classdef BloodVessel < PeriodicCylinder & ForceMixin
    methods
        function obj = BloodVessel(n, m, varargin)
            obj@PeriodicCylinder(n, m);
            obj@ForceMixin(varargin{:});
        end

        function x = shape(~, params)
            theta = params(:, 1);
            zeta = params(:, 2) / (2 * pi);

            rho = @(x) tanh(10*(x+0.5))-tanh(10*(x-0.5));
            d = 0;
            for n = -1:1
                d = d + rho(-1 + 2 * (zeta - n));
            end

            x = [2^-9 + 15e-4 * (1 - 1/12 * d) .* cos(theta), ...
                 2^-8 * zeta, ...
                 2^-9 + 15e-4 * (1 - 1/12 * d) .* sin(theta)];
        end
    end

%    methods(Static)
%        function params = sample(n)
%            rho = @(x) tanh(10*(x+0.5))-tanh(10*(x-0.5));
%            drho = @(x) -10 * tanh(10*(x+0.5)).^2 + 10 * tanh(10*(x-0.5)).^2;
%            r = 0;
%            dr = 0;
%
%            sc = 2 * pi / 6;
%            layers = floor(sqrt(n / sc));
%            uniform = ((1:n)-0.5)' / n;
%            for n = -1:1
%                r = r + rho(-1 + 2 * (uniform - n));
%                dr = dr + drho(-1 + 2 * (uniform - n));
%            end
%            % This spacing isn't quite right (constants should take into
%            % account the different radii).
%            spacing = (1 + dr.^2).^(-1/4);
%            scaling = 1 / sum(spacing);
%            z = cumsum(scaling * spacing) - scaling * spacing(1);
%            theta = mod(2 * pi * layers * uniform, 2 * pi);
%            params = [theta z];
%        end
%    end
end
