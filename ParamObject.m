classdef ParamObject < handle
    properties(Access = protected)
        params
    end

    methods
        function obj = ParamObject(varargin)
            obj.params = varargin;
        end
    end
end
