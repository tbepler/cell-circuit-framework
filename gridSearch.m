function V = gridSearch( fnc, varargin )
    if length( varargin ) >= 1
        V = zeros( size( varargin{1} ) );
        %args = zeros( length( varargin ), 1 );
        parfor i = 1:numel(varargin{1})
            args = zeros( length( varargin ), 1 );
            for j = 1:length(args)
                args(j) = varargin{j}(i);
            end
            V(i) = fnc( args );
        end
            
    else
        V = fnc();
    end
        
        
end