classdef Gene < handle
    
    properties(SetAccess = immutable)
        name
        mRNA
    end
    
    properties(SetObservable = true)
        kOnBasal
        kOffBasal
    end
    
    methods
        function obj = Gene( name, mRNA, kOnBasal, kOffBasal )
            obj.name = name;
            obj.mRNA = mRNA;
            obj.kOnBasal = kOnBasal;
            obj.kOffBasal = kOffBasal;
        end
        
        function self = accept( self, sys, initVal )
            %add compositor for the protein
            comp = sys.AddCompositor( self.name, initVal );
            
            %add rate constants
            sys.AddConstant( kOnBasalName(self), self.kOnBasal );
            sys.AddConstant( kOffBasalName(self), self.kOffBasal );
            
            %add basal transcription model
            rnaComp = sys.getCompositor( self.mRNA );
            rnapComp = sys.RNAP;
            
            complex = sys.AddCompositor( rnapComplexName(self), 0 );
            binding = [ kOnBasalName(self), ' * ', self.name ...
                , ' * ', sys.RNAP_NAME ];
            unbinding = [ kOffBasalName(self) '*' ...
                rnapComplexName(self) ];
            txn = [ sys.K_TXN_NAME '*' rnapComplexName(self) ];
            
            sys.AddPart( Part( [self.name ' transcription'] ...
                , [ comp rnapComp complex rnaComp ] ...
                , [ Rate( [ '-' binding '+' unbinding '+' txn ] ) ...
                    Rate( [ '-' binding '+' unbinding '+' txn ] ) ...
                    Rate( [ '+' binding '-' unbinding '-' txn ] ) ...
                    Rate( [ '+' txn ] ) ] ) );
            
        end
        
        function self = update( self, sys )
            sys.ChangeConstant( kOnBasalName(self), self.kOnBasal );
            sys.ChangeConstant( kOffBasalName(self), self.kOffBasal );
        end
        
    end 
    
end

function name = rnapComplexName( obj )
    name = ['Ribo',obj.name];
end

function name = kOnBasalName( obj )
    name = [ 'k_onBasal_',obj.name ];
end

function name = kOffBasalName( obj )
    name = [ 'k_offBasal_',obj.name ];
end