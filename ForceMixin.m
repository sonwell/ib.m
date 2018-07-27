classdef (HandleCompatible) ForceMixin
    methods(Abstract)
        shape(self, params)
    end

    properties
        force_handle
    end

    methods
        function self = ForceMixin(varargin)
            handle = CombineForces(varargin{:});
            self.force_handle = handle(self);
        end

        function f = force(self)
            f = self.force_handle();
        end
    end
end
