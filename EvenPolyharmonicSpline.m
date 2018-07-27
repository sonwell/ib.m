classdef EvenPolyharmonicSpline
    properties
        n
    end

    methods
        function obj = EvenPolyharmonicSpline(n)
            obj.n = n;
        end

        function m = phi(obj, r)
            meps = 2^-52;
            n = obj.n;

            m = r.^n .* log(r + meps);
        end

        function m = dphi(obj, r, rdr)
            meps = 2^-52;
            n = obj.n;

            m = r.^(n-2) .* (1 + n * log(r + meps)) .* rdr;
        end

        function m = ddphi(obj, r, drdr, rdr2)
            meps = 2^-52;
            n = obj.n;

            m = obj.dphi(r, drdr) + r.^(n-4) .* (2 * n - 2 + n * (n-2) * log(r + meps)) .* rdr2;
        end
    end
end
