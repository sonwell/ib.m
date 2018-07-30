function result = indexfn(fn, dims)
    arity = nargin(fn);
    output = cell(arity, 1);
    [output{:}] = ndgrid((1:dims)');
    result = arrayfun(fn, output{:}, 'UniformOutput', false);
end
