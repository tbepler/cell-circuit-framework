classdef mRNA < handle
    
    properties(SetAccess = immutable)
        name
        protein
    end
    
    properties(SetObservable = true)
        kDeg
        kOnRibo
        kOffRibo
    end
    
    methods
        function obj = mRNA( name, protein, kDeg, kOnRibo, kOffRibo )
            obj.name = name;
            obj.protein = protein;
            obj.kDeg = kDeg;
            obj.kOnRibo = kOnRibo;
            obj.kOffRibo = kOffRibo;
        end
        
        function self = accept( self, sys, initVal )
            %add compositor for the protein
            comp = sys.AddCompositor( self.name, initVal );
            
            %add rate constants
            sys.AddConstant( kDegName(self) , self.kDeg );
            sys.AddConstant( kOnRiboName(self), self.kOnRibo );
            sys.AddConstant( kOffRiboName(self), self.kOffRibo );
            
            %add part for degradation
            rate = Rate( [ '- ', kDegName(self), ' * ', self.name ] );
            sys.AddPart( Part( [self.name, ' degradation'] ...
                , comp, rate ) );
            
            %add model for translation
            protComp = sys.getCompositor( self.protein );
            riboComp = sys.Ribo;
            
            complex = sys.AddCompositor( riboComplexName(self), 0 );
            binding = [ kOnRiboName(self), ' * ', self.name ...
                , ' * ', sys.RIBO_NAME ];
            unbinding = [ kOffRiboName(self) '*' ...
                riboComplexName(self) ];
            tln = [ sys.K_TLN_NAME '*' riboComplexName(self) ];
            
            sys.AddPart( Part( [self.name ' translation'] ...
                , [ comp riboComp complex protComp ] ...
                , [ Rate( [ '-' binding '+' unbinding '+' tln ] ) ...
                    Rate( [ '-' binding '+' unbinding '+' tln ] ) ...
                    Rate( [ '+' binding '-' unbinding '-' tln ] ) ...
                    Rate( [ '+' tln ] ) ] ) );
            
        end
        
        function self = update( self, sys )
            sys.ChangeConstant( kDegName(self), self.kDeg );
            sys.ChangeConstant( kOnRiboName(self), self.kOnRibo );
            sys.ChangeConstant( kOffRiboName(self), self.kOffRibo );
        end
        
    end 
    
end

function name = riboComplexName( obj )
    name = ['Ribo',obj.name];
end

function name = kDegName( obj )
    name = [ 'gamma_',obj.name ];
end

function name = kOnRiboName( obj )
    name = [ 'k_onRibo_',obj.name ];
end

function name = kOffRiboName( obj )
    name = [ 'k_offRibo_',obj.name ];
end