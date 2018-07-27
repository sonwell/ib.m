function fn = CombineForces(varargin)
    if isempty(varargin)
        fn = @null_forces;
    elseif length(varargin) == 1
        fn = varargin{1};
    else
        fn = combine(varargin{1}, CombineForces(varargin{2:end}));
    end
end

function fn = null_forces(obj)
    gr = obj.geometry('reference');
    grs = gr.at_sample_sites;
    f0 = 0 * [grs.r];
    function f = zero()
        f = f0;
    end
    fn = @zero;
end

function fn = combine(left, right)
    function fn = combined_initializer(obj)
        left_force = left(obj);
        right_force = right(obj);

        function f = forces()
            f = left_force() + right_force();
        end
        fn = @forces;
    end
    fn = @combined_initializer;
end
