function fn = CombineForces(varargin)
    if length(varargin) == 0
        fn = @null_forces;
    elseif length(varargin) == 1
        fn = varargin{1};
    else
        fn = combine(varargin{1}, CombineForces(varargin{2:end}));
    end
end

function fn = null_forces(orig)
    function f = zero(curr)
        f = 0 * curr.E;
    end
end

function fn = combine(left, right)
    function fn = combined_initializer(orig)
        left_force = left(orig);
        right_force = right(orig);

        function f = forces(curr)
            f = left_force(curr) + right_force(curr);
        end
    end
    fn = @combined_initializer;
end
