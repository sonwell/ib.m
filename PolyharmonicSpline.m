function phs = PolyharmonicSpline(n)
    if mod(n, 2) == 0
        phs = EvenPolyharmonicSpline(n);
    else
        phs = OddPolyharmonicSpline(n);
    end
end
