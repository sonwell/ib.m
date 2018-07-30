classdef Operators
    properties
        id
        d
        dd
        int
    end

    methods
        function obj = Operators(data, sample, itp, qdr, metric, rbf, poly, weights, qpoly)
            n = size(data, 1);
            dims = size(data, 2);
            trim = @(m) m(:, 1:n);

            [r, rdr, drdr] = metric(data, sample);
            psi = rbf.phi(r);
            dpsi = indexfn(@(i) rbf.dphi(r, rdr{i}), dims);
            ddpsi = indexfn(@(i, j) rbf.ddphi(r, drdr{i, j}, rdr{i} .* rdr{j}), dims);
            q = poly.p(sample);
            dq = indexfn(@(i) poly.dp(sample, i), dims);
            ddq = indexfn(@(i, j) poly.dp(sample, [i j]), dims);

            obj.id = trim((itp' \ [psi q]')');
            obj.d = indexfn(@(i) trim((itp' \ [dpsi{i} dq{i}]')'), dims);
            obj.dd = indexfn(@(i, j) trim((itp' \ [ddpsi{i, j} ddq{i, j}]')'), dims);

            id2 = trim((qdr' \ [psi qpoly.p(sample)]')');
            sc = sum(id2, 1);
            sp = bsxfun(@rdivide, id2, sc);
            obj.int = weights * sp';
        end
    end
end
