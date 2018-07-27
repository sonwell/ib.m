classdef ClosedSurface < ParamObject
    % ClosedSurface Geometric representation of a closed shape in 3D.

    properties (Dependent)
        x
    end

    properties (SetAccess = private)
        data_sites  % Interpolation point parameters
        sample_sites  % Evaluation point parameters
        rbf
        poly
        qpoly

        phi
        p
        q

        opd
        ops
    end

    properties (Access = protected)
        gd
        gs

        gdr
        gsr
    end

    methods (Abstract)
        shape(self, params)
    end

    methods (Abstract, Static)
        sample(n)
        metric(data, sample)
    end

    methods
        function self = ClosedSurface(n, m, rbf, poly, qpoly, rhs, wfn, varargin)
            self@ParamObject(varargin{:});
            data = self.sample(n);
            sample = self.sample(m);
            self.data_sites = data;
            self.sample_sites = sample;

            self.rbf = rbf;
            self.poly = poly;
            self.qpoly = qpoly;

            r = self.metric(data, data);
            self.phi = self.rbf.phi(r);
            self.p = self.poly.p(data);
            self.q = self.qpoly.p(data);
            p0 = zeros(size(self.p, 2));
            q0 = zeros(size(self.q, 2));

            itp = [self.phi self.p; self.p' p0];
            qdr = [self.phi self.q; self.q' q0];
            w = wfn(data);
            ds = self.trim((qdr' \ [zeros(1, n) rhs(:)']')') ./ w(:)';

            self.opd = Operators(data, data, itp, qdr, @self.metric, rbf, poly, ds, qpoly);
            self.ops = Operators(data, sample, itp, qdr, @self.metric, rbf, poly, ds, qpoly);
            self.opd.id = 1;

            self = self.update_geometry(self.shape(data));
            self.gdr = self.gd;
            self.gsr = self.gs;
        end

        function m = trim(self, M)
            n = size(self.data_sites, 1);
            m = M(:, 1:n);
        end

        function x = get.x(self)
            x = self.gd.r;
        end

        function self = set.x(self, x)
            self = self.update_geometry(x);
        end

        function [px, py, pz] = surf(self)
            data = self.data_sites;
            fill = zeros(size(self.p, 2));
            itp = [self.phi self.p; self.p' fill];
            w = itp \ [self.x; zeros(size(self.p, 2), size(self.x, 2))];

            function y = surfer(u, v, c)
                sample = [u(:) v(:)];
                r = self.metric(data, sample);
                psi = self.rbf.phi(r);
                s = self.poly.p(sample);
                y = [psi s] * c;
            end
            
            px = @(u, v) reshape(surfer(u, v, w(:, 1)), size(u));
            py = @(u, v) reshape(surfer(u, v, w(:, 2)), size(u));
            pz = @(u, v) reshape(surfer(u, v, w(:, 3)), size(u));
        end

        function g = geometry(self, which)
            g = feval(which, self);
        end
    end

    methods (Access = private)
        function self = update_geometry(self, x)
            self.gd = Geometry(x, self.opd);
            self.gs = Geometry(x, self.ops);
        end

        function g = reference(self)
            g.at_sample_sites = self.gsr;
            g.at_data_sites = self.gdr;
        end

        function g = current(self)
            g.at_sample_sites = self.gs;
            g.at_data_sites = self.gd;
        end
    end
end
