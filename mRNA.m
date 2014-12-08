classdef mRNA < handle
    
    properties(SetAccess = immutable)
        name
        protein
    end
    
    properties(SetObservable = true)
        kDeg
        kOnRibo
        kOffRibo
        kRiboInit
    end
    
    methods
        function obj = mRNA( name, protein, kDeg, kOnRibo, kOffRibo, kRiboInit )
            %args are ( name, protein, degredation rate, ribosome on rate,
            %ribosome off rate )
            obj.name = name;
            obj.protein = protein;
            obj.kDeg = kDeg;
            obj.kOnRibo = kOnRibo;
            obj.kOffRibo = kOffRibo;
            obj.kRiboInit = kRiboInit;
        end
        
        function self = accept( self, sys, initVal )
            %add compositor for the protein
            comp = sys.AddCompositor( self.name, initVal );
            
            %add rate constants
            sys.AddConstant( kDegName(self) , self.kDeg );
            sys.AddConstant( kOnRiboName(self), self.kOnRibo );
            sys.AddConstant( kOffRiboName(self), self.kOffRibo );
            sys.AddConstant( kRiboInitName(self), self.kRiboInit );
            
            %add part for degradation
            deg = Rate( [ '- ', kDegName(self), ' * ', self.name ] );
            sys.AddPart( Part( [self.name, ' degradation'] ...
                , comp, deg ) );
            
            %add model for translation
            protComp = sys.getCompositor( self.protein );
            riboComp = sys.Ribo;
            
            complex = sys.AddCompositor( riboComplexName(self), 0 );
            binding = [ kOnRiboName(self), ' * ', self.name ...
                , ' * ', sys.RIBO_NAME ];
            unbinding = [ kOffRiboName(self) '*' ...
                riboComplexName(self) ];
            tln = [ kRiboInitName(self) '*' riboComplexName(self) ];
            
            sys.AddPart( Part( [self.name ' translation'] ...
                , [ comp riboComp complex protComp ] ...
                , [ Rate( [ '-' binding '+' unbinding '+' tln ] ) ...
                    Rate( [ '-' binding '+' unbinding '+' tln ] ) ...
                    Rate( [ '+' binding '-' unbinding '-' tln ] ) ...
                    Rate( [ '+' tln ] ) ] ) );
                
                %Ribosome bound mRNA also degrades
             deg = ['- ', kDegName(self), ' * ', riboComplexName(self) ];
             sys.AddPart( Part( [riboComplexName(self), ' mRNA degredation'] ...
                 , [ complex riboComp ] ...
                 , [ Rate( deg ) ...
                     Rate( [ '-' deg ] ) ] ) );
            
        end
        
        function self = update( self, sys )
            sys.ChangeConstant( kDegName(self), self.kDeg );
            sys.ChangeConstant( kOnRiboName(self), self.kOnRibo );
            sys.ChangeConstant( kOffRiboName(self), self.kOffRibo );
            sys.ChangeConstant( kRiboInitName(self), self.kRiboInit );
        end
        
    end 
    
end

function name = riboComplexName( obj )
    name = ['Ribo',obj.name];
end

function name = kDegName( obj )
    name = [ 'gamma_',obj.name ];
end

function name = kRiboInitName( obj )
    name = [ 'k_riboInit_', obj.name ];
end

function name = kOnRiboName( obj )
    name = [ 'k_onRibo_',obj.name ];
end

function name = kOffRiboName( obj )
    name = [ 'k_offRibo_',obj.name ];
end