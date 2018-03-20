classdef ForceMixin
    properties
        force_handle
    end

    methods
        function obj = ForceMixin(params, varargin)
            handle = CombineForces(varargin{:});
            x = obj.shape(params);
            orig = obj.geometry(x);
            obj.force_handle = handle(orig);
        end

        function [f, y] = force(obj, x)
            curr = obj.geometry(x);
            f = obj.force_handle(curr);
            y = curr.x;
        end
    end
end