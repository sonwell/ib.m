classdef ForceMixin
    methods(Abstract, Static)
        shape(params)
    end

    methods(Abstract)
        geometry(obj, x)
    end

    properties(Abstract)
        data_sites
    end

    properties
        force_handle
    end

    methods
        function obj = ForceMixin(varargin)
            handle = CombineForces(varargin{:});
            x = obj.shape(obj.data_sites);
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
