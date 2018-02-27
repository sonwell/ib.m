function calculator = ForceCalculator(structure)
    handle = structure.force;
    transform = structure.reference;
    params = structure.interpolation.parametrization;
    id = structure.id;
    
    X0 = transform(params);
    geometry = SurfaceGeometryCalculator(structure);
    orig = geometry(X0);
    force = handle(orig);

    function [f, x] = f_calculator(X)
        curr = geometry(X);
        f = force(curr);
        x = id * X;
    end
    calculator = @f_calculator;
end
