classdef OddPolyharmonicSpline
    properties
        n
    end

    methods
        function obj = OddPolyharmonicSpline(n)
            obj.n = n;
        end

        function m = phi(obj, r)
            n = obj.n;

            m = r.^n;
        end

        function m = dphi(obj, r, rdr)
            n = obj.n;

            m = n * r.^(n-2) .* rdr;
        end

        function m = ddphi(obj, r, drdr, rdr2)
            n = obj.n;

            m = obj.dphi(r, drdr) + n * (n-2) * r.^(n-4) .* rdr2;
        end
    end
end
