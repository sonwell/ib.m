classdef Geometry
    properties
        r
        dr
        ddr
        n

        E
        F
        G

        e
        f
        g

        H
        K

        I
        II

        ds
    end

    methods
        function self = Geometry(x, ops)
            dims = 2;

            self.r = ops.id * x;
            self.dr = indexfn(@(i) ops.d{i} * x, dims);
            self.ddr = indexfn(@(i, j) ops.dd{i, j} * x, dims);

            self.E = dot(self.dr{1}, self.dr{1}, 2);
            self.F = dot(self.dr{1}, self.dr{2}, 2);
            self.G = dot(self.dr{2}, self.dr{2}, 2);
            self.I = self.E .* self.G - self.F .^ 2;
            detg = sqrt(self.I);

            self.n = bsxfun(@rdivide, cross(self.dr{1}, self.dr{2}, 2), detg);

            self.e = dot(self.ddr{1, 1}, self.n, 2);
            self.f = dot(self.ddr{1, 2}, self.n, 2);
            self.g = dot(self.ddr{2, 2}, self.n, 2);
            self.II = self.e .* self.g - self.f .^ 2;

            self.H = (self.e .* self.G - 2 * self.f .* self.F + self.E .* self.g) ./ (2 * self.I);
            self.K = self.II ./ self.I;

            self.ds = ops.int .* detg';
        end
    end
end
