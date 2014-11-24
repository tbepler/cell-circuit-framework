function [xelems, yelems] = listDiff( x, y )

    xelems = x;
    for elem = y
        i = firstIndexOf( xelems, elem );
        if i > 0
            xelems(i) = [];
        end
    end
    
    yelems = y;
    for elem = x
        i = firstIndexOf( yelems, elem );
        if i > 0
            yelems(i) = [];
        end
    end

%     xmap = counter( x );
%     ymap = counter( y );
%     
%     for k = keys( xmap )
%         xcount = lookup( xmap, k, 0 );
%         if xcount > 0
%             ycount = lookup( ymap, k, 0 );
%             if ycount > 0
%                 remove = min( [xcount, ycount] );
%                 xmap(k) = xcount - remove;
%                 ymap(k) = ycount - remove;
%             end
%         end
%     end
%     
%     xelems = asVector( xmap );
%     yelems = asVector( ymap );

end

function idx = firstIndexOf( v, e )
    for idx = 1 : numel( v )
        if v(idx ) == e
            return
        end
    end
    idx = 0;
end

function v = asVector( counter )
    %len = sumCounter( counter );
    v = [];
    i = 1;
    for k = keys( counter )
        n = counter( k );
        for j = 1 : n
            v(i) = k;
            i = i + 1;
        end
    end
end

function s = sumCounter( map )
    s = 0;
    for k = keys( map )
        count = map( k );
        s = s + count;
    end
end

function val = lookup( map, k, default )
    if isKey( map, k )
        val = map(k);
    else
        val = default;
    end
end