classdef PolynomialSubset
    properties
        projection
        dimensions
        degree
        exponents
    end

    methods
        function obj = PolynomialSubset(degree, dimensions, keep)
            trim = @(m) m(:, 2:end);
            remove = ones(dimensions, 1);
            remove(keep(:)) = 0;
            I = speye(dimensions);
            projection = I(remove == 0, :);

            obj.projection = projection;
            obj.dimensions = dimensions;
            obj.degree = degree;
            obj.exponents = trim(obj.multiindex(size(projection, 1) + 1, degree));
        end

        function P = p(obj, x)
            e = obj.exponents;
            y = x * obj.projection';
            nx = size(y, 1);
            ne = size(e, 1);
            [xi, ei] = ndgrid(1:nx, 1:ne);

            P = reshape(prod(y(xi(:), :) .^ e(ei(:), :), 2), nx, ne);
        end

        function P = dp(obj, x, d)
            e = obj.exponents;
            y = x * obj.projection';
            a = zeros(1, obj.dimensions);
            acc = accumarray(d(:), 1);
            a(1:numel(acc)) = acc;
            c = obj.coefficients(d);
            n = ones(size(e, 1), 1) * (a * obj.projection');
            e = max(0, e - n);
            nx = size(y, 1);
            ne = size(e, 1);
            [xi, ei] = ndgrid(1:nx, 1:ne);

            P = reshape(c(ei(:)) .* prod(y(xi(:), :) .^ e(ei(:), :), 2), nx, ne);
        end
    end

    methods(Access = private)
        function c = coefficients(obj, d)
            p = obj.exponents * obj.projection;
            c = ones(size(p, 1), 1);
            for j = d
                c = c .* p(:, j);
                p(:, j) = p(:, j) - 1;
            end
        end

        function r = multiindex(obj, d, l)
            if d == 1
                r = l;
                return;
            end

            offset = 0;
            r = zeros(nchoosek(d-1+l, l), d);
            for n=l:-1:0
                s = obj.multiindex(d-1, l-n);
                ss = size(s, 1);
                o = n * ones(ss, 1);
                r((offset+1):(offset+ss), :) = [o s];
                offset = offset + ss;
            end
        end
    end
end
