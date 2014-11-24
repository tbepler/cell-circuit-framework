classdef Protein < handle
    
    properties(SetAccess = immutable)
        name
    end
    
    properties(SetObservable = true)
        kDeg
    end
    
    methods
        function obj = Protein( name, kDeg )
            %args are (name, degredation rate)
            obj.name = name;
            obj.kDeg = kDeg;
        end
        
        function self = accept( self, sys, initVal )
            %add compositor for the protein
            compP = sys.AddCompositor( self.name, initVal );
            
            %add degredation rate constant
            sys.AddConstant( kDegName(self) , self.kDeg );
            
            %add part for degradation
            rate = [ '- ', kDegName(self), ' * ', self.name ];
            sys.AddPart( Part( [self.name, ' degradation'] ...
                , compP, Rate( rate ) ) );
        end
        
        function self = update( self, sys )
            sys.ChangeConstant( kDegName(self), self.kDeg );
        end
        
    end 
end

function name = kDegName( p )
    name = [ 'gamma_',p.name ];
end