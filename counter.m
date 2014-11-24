function map = counter( xs )
    map = containers.Map();
    for x = xs
        if isKey( map, x )
            map(x) = map(x) + 1;
        else
            map(x) = 1;
        end
    end
end
    