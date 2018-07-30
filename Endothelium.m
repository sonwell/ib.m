classdef Endothelium < PeriodicSheet & ForceMixin
    methods
        function obj = Endothelium(n, m, varargin)
            obj@PeriodicSheet(n, m);
            obj@ForceMixin(varargin{:});
        end

        function x = shape(~, params)
            a = params(:, 1) / (2 * pi);
            b = params(:, 2) / (2 * pi);

            x = [2^-8*a, 2^-8*b, 0*a];
        end
    end
end
