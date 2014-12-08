classdef Reactions
    %REACTIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods (Static)
        
        function fnc = complexAtEq( X_tot, Y_tot, K_D )
            syms complex_temp;
            slns = solve( X_tot * ( Y_tot - complex_temp ) / ( K_D + Y_tot - complex_temp ) ...
                == complex_temp, complex_temp );
%             slns = solve( k_on * ( X_tot - complex_temp ) * ( Y_tot - complex_temp ) ...
%                 == k_off * complex_temp, complex_temp );
            fnc = slns(2);
        end
        
        function fnc = reactionComplexAtEq( X_tot, Y_tot, K_M )
            syms complex_temp;
            slns = solve( X_tot * ( Y_tot - complex_temp ) / ( K_M + Y_tot - complex_temp ) ...
                == complex_temp, complex_temp );
            fnc = slns(2);
        end
        
        function fnc = mrnaProduction( promoter, RNAP, K_M, k_ocl )
            complex = Reactions.complexAtEq( promoter, RNAP, K_M );
            fnc = k_ocl * complex;
        end
        
        function fnc = mrnaAtEq( production_terms, k_deg )
            P = 0;
            for i = 1 : length(production_terms)
                P = P + production_terms(i);
            end
            fnc = P / k_deg;
        end
        
    end
    
end

