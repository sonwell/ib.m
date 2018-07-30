function RBFFluidTest(domain, rho, mu, cell, vessel)
    nx = domain.nx;
    ny = domain.ny;
    nz = domain.nz;
    h = (domain.bounds(:, 2) - domain.bounds(:, 1))' ./ [nx ny nz];

    volume = prod(h);

    Ue = zeros(nx * ny * nz, 3);
    Ff = zeros(nx * ny * nz, 3);

    offsets = [0.0, 0.5, 0.5; 0.5, 0.0, 0.5; 0.5, 0.5, 0.0];  % MAC offsets
    % Helpers functions for spreading and interpolation
    spread = @(x, f, da) ...
        [Transfer('l2e', domain, x, offsets(1, :), f(:, 1) .* da), ...
         Transfer('l2e', domain, x, offsets(2, :), f(:, 2) .* da), ...
         Transfer('l2e', domain, x, offsets(3, :), f(:, 3) .* da)];
    interpolate = @(x, u, dv) ...
        [Transfer('e2l', domain, x, offsets(1, :), u(:, 1) .* dv), ...
         Transfer('e2l', domain, x, offsets(2, :), u(:, 2) .* dv), ...
         Transfer('e2l', domain, x, offsets(3, :), u(:, 3) .* dv)];

    k = 1e-6; %0.015 * exp(mean(log(h)));
    [solver, ~, handles] = NavierStokesSolver(domain, rho, mu, k);

    Uy = 0;
    Ue(:, 2) = Uy;
    Fy = 8000; %-mu * (L * Ue(:, 2));
    Ff(:, 2) = Fy;

    avg = @(X, w) (w * X) ./ sum(w);

    grc = cell.geometry('reference');
    grv = vessel.geometry('reference');

    dsc = grc.at_sample_sites.ds;
    dsv = grv.at_sample_sites.ds;

    force_info = @(obj) deal(obj.force(), obj.ops.id * obj.x);

    fig = figure;
    for s = 0:150000
        xcd = cell.x;
        xvd = vessel.x;

        % Backward Euler step
        [fc, xcs] = force_info(cell);
        [fv, xvs] = force_info(vessel);
        Fc = spread(xcs, fc, dsc');
        Fv = spread(xvs, fv, 1);
        [Us, ~] = solver(handles{1}, Ue, Ue, Fc + Fv + Ff);
        uc = interpolate(xcd, Us, volume);
        uv = interpolate(xvd, Us, volume);
        cell.x = xcd + (k/2) * uc;
        vessel.x = xvd + (k/2) * uv;

        % Crank-Nicolson step
        [fc, xcs] = force_info(cell);
        [fv, xvs] = force_info(vessel);
        Fc = spread(xcs, fc, dsc');
        Fv = spread(xvs, fv, 1);
        [Ue, ~] = solver(handles{2}, Ue, Us, Fc + Fv + Ff);
        uc = interpolate(xcd, Ue, volume);
        uv = interpolate(xvd, Ue, volume);
        cell.x = xcd + k * uc;
        vessel.x = xvd + k * uv;

        fprintf('%d: %12g %12g %12g %12g %12g %12g\n', s, avg(xcs, dsc), avg(fv, dsv));
        if mod(s, 100) == 0
            dlmwrite(sprintf('data/rbc_%04d.txt', s/100), cell.x, 'delimiter', '\t', 'precision', 18);
            dlmwrite(sprintf('data/vessel_%04d.txt', s/100), vessel.x, 'delimiter', '\t', 'precision', 18);
        end
        if mod(s, 1000) == 0
            plotter(domain, s, k, cell, vessel);

            set(fig, 'Units', 'Inches');
            pos = get(fig, 'Position');
            set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
            saveas(fig, sprintf('simulation%03d.pdf', s/10));
        end
    end
    close all;
end

function plotter(domain, s, k, cell, vessel)
    lower = domain.bounds(:, 1)';
    length = domain.bounds(:, 2)' - lower;

    function psurf(x, y, z, n, params, fc)
        m = n(1, :);
        nr = n(2, :) - m;
        [dx, dy, dz] = ndgrid(0:nr(1), 0:nr(2), 0:nr(3));

        for i = 1:numel(dx)
            sx = @(u, v) x(u, v) - (m(1) - dx(i)) * length(1);
            sy = @(u, v) y(u, v) - (m(2) - dy(i)) * length(2);
            sz = @(u, v) z(u, v) - (m(3) - dz(i)) * length(3);
            srf = ezsurf(sx, sy, sz, params);
            srf.FaceColor = fc;
            srf.EdgeColor = 'none';
        end
    end

    function n = frames(x)
        fr = floor(bsxfun(@rdivide, bsxfun(@minus, x, lower), length));
        n = [min(fr); max(fr)];
    end

    function subplotter(vp, cp)
        hold on;
        psurf(vx, vy, vz, vn, vp, 'b');
        psurf(cx, cy, cz, cn, cp, 'r');
        hold off;
        axis equal;
        xlim(domain.bounds(1, :));
        ylim(domain.bounds(2, :));
        zlim(domain.bounds(3, :));
        title(sprintf('$t = %f$', s * k), 'Interpreter', 'latex');
        camlight;
    end

    [cx, cy, cz] = cell.surf();
    [vx, vy, vz] = vessel.surf();
    cn = frames(cell.x);
    vn = frames(vessel.x);

    % Plot
    %subplot(1, 2, 1);
    subplotter([0, 2*pi, 0, 2*pi], [0, pi, 0, 2*pi]);
    view([-90 0]);

    %subplot(1, 2, 2);
    %subplotter([-pi, pi, 0, 2*pi], [0, pi, 0, 2*pi]);
    %view([0 0]);
end
