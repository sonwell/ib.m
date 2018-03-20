classdef PolyharmonicSpline
    properties
        n
    end

    methods
        function obj = PolyharmonicSpline(n)
            obj.n = n;
        end

        function m = phi(obj, r)
            meps = 2^-52;
            n = obj.n;

            m = r.^n .* log(r + meps);
        end

        function m = dphi(obj, r, dr)
            meps = 2^-52;
            n = obj.n;

            m = r.^(n-2) .* (1 + n * log(r + meps)) .* dr;
        end

        function m = ddphi(obj, r, d2r, dr2)
            meps = 2^-52;
            n = obj.n;

            m = obj.dphi(r, d2r) + r.^(n-4) .* (2 * n - 2 + n * (n-2) * log(r + meps)) .* dr2;
        end
    end
end
