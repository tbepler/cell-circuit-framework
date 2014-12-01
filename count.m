function [elems, counts] = count( X )
    
    elems = unique( X );
    counts = zeros( length(elems), 1 );
    for i = 1:length(elems)
        counts(i) = sum( X == elems(i) );
    end

end