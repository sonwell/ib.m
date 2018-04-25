classdef Polynomials < PolynomialSubset
    methods
        function obj = Polynomials(degree, dimensions)
            obj@PolynomialSubset(degree, dimensions, 1:dimensions);
        end
    end
end
