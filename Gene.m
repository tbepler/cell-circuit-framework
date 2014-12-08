classdef Gene < handle
    
    properties(SetAccess = immutable)
        name
        mRNA
        promoterStates
        rnapOn
        rnapOff
        rnapInit
        transitionMatrix
    end
    
    properties(SetObservable = true)
        %kOnBasal
        %kOffBasal
    end
    
    methods
        %promoterStates is cell array where each entry is a set of proteins
        %FIRST STATE IS NAKED PROMOTER STATE
        %rnapOn is vector of same length as promoterStates specifying the
        %on rate of RNAP
        %rnapOff "" off rate of RNAP
        %transition matrix is an NxN matrix where N is the number of
        %promoter states specifying the transition rate for going from
        %state i to state j
        function obj = Gene( name, mRNA, promoterStates, rnapOn, rnapOff, rnapInit, transitionMatrix )
            obj.name = name;
            obj.mRNA = mRNA;
            obj.promoterStates = promoterStates;
            obj.rnapOn = rnapOn;
            obj.rnapOff = rnapOff;
            obj.rnapInit = rnapInit;
            obj.transitionMatrix = transitionMatrix;
            %obj.kOnBasal = kOnBasal;
            %obj.kOffBasal = kOffBasal;
        end
        
        function self = accept( self, sys, initVal )
            %add compositor for the protein
            %comp = sys.AddCompositor( self.name, initVal );
            
            %add promoter state compositors to the system
            N = length( self.promoterStates );
            stateCompositors = cell( N, 1 );
            stateNames = cell( N, 1 );
            for i = 1 : N
                proteins = self.promoterStates{i};
                sName = self.name;
                for p = proteins
                    sName = [ sName, p.name ];
                end
                stateNames{ i } = sName;
                %define the RNAP unbound and bound compositors for this
                %state
                if i == 1
                    sComp = sys.AddCompositor( sName, initVal );
                else
                    sComp = sys.AddCompositor( sName, 0 );
                end
                stateCompositors{ i } = sComp;
                sRNAPName = [ sName, sys.RNAP_NAME ];
                sRNAPComp = sys.AddCompositor( sRNAPName, 0 );
                
                %define RNAP on and off rates for this state
                rnapOnConst = [ 'k_', sys.RNAP_NAME, '_on_', sName ];
                sys.AddConstant( rnapOnConst, self.rnapOn(i) );
                rnapOffConst = [ 'k_', sys.RNAP_NAME, '_off_', sRNAPName ];
                sys.AddConstant( rnapOffConst, self.rnapOff(i) );
                rnapInitConst = [ 'k_', sys.RNAP_NAME, '_init_', sRNAPName ];
                sys.AddConstant( rnapInitConst, self.rnapInit(i) );
                
                %add parts for RNAP binding and transcription
                rna = sys.getCompositor( self.mRNA );
                rnap = sys.RNAP;
                binding = [ rnapOnConst, ' * ', sName, ' * ', sys.RNAP_NAME];
                unbinding = [ rnapOffConst, ' * ', sRNAPName ];
                txn = [ rnapInitConst, ' * ', sRNAPName ];
                
                sys.AddPart( Part( [sName, ' transcription'] ...
                    , [ sComp, rnap, sRNAPComp, rna ] ...
                    , [ Rate( [ '-' binding '+' unbinding '+' txn ] ) ...
                        Rate( [ '-' binding '+' unbinding '+' txn ] ) ...
                        Rate( [ '+' binding '-' unbinding '-' txn ] ) ...
                        Rate( [ '+' txn ] ) ] ) );
                
            end
                
            %add transition parts for transitioning between promoter states
            for i = 1 : N
                fromName = stateNames{ i };
                fromComp = stateCompositors{ i };
                fromProts = self.promoterStates{ i };
                for j = 1 : N
                    if i ~= j
                        toName = stateNames{ j };
                        toComp = stateCompositors{ j };
                        toProts = self.promoterStates{ j };
                        rate = [ 'k_', fromName, '_', toName ];
                        sys.AddConstant( rate, self.transitionMatrix( i, j ) );
                        
                        [removeProts, addProts] = listDiff( fromProts, toProts );
                        [removeProts, removeCoops] = count( removeProts );
                        [addProts, addCoops] = count( addProts );
                        trans = [ 'exp(' rate, '+ log(', fromName ')' ];
                        A = length( addProts );
                        R = length( removeProts);
                        for k = 1:A
                            trans = [ trans, '+ log(', addProts(k).name, ')*', num2str(addCoops(k)) ];
                        end
                        trans = [ trans ')' ];
                        %disp( trans );
                        compositors = [ fromComp, toComp ];
                        rates = [ Rate( [ '-' trans ] ) ...
                                 ,Rate( [ '+' trans ] ) ];
                        %add the removed protein compositors and rates
                        for k = 1:R
                            compositor = sys.getCompositor( removeProts(k) );
                            compositors = [ compositors, compositor ];
                            rates = [ rates ...
                                , Rate( [ '+' num2str(removeCoops(k)) '*' trans ] ) ];
                        end
                        %add the added protein compositors and rates
                        for k = 1:A
                            compositor = sys.getCompositor( addProts(k) );
                            compositors = [ compositors, compositor ];
                            rates = [ rates ...
                                , Rate( [ '-' num2str(addCoops(k)) '*' trans ] ) ] ;
                        end
                        
                        %add the part
                        sys.AddPart( Part( [ fromName, '->', toName ] ...
                            , compositors ...
                            , rates ) );
                        
                    end
                end
            end
            
        end
        
        function self = update( self, sys )
            %sys.ChangeConstant( kOnBasalName(self), self.kOnBasal );
            %sys.ChangeConstant( kOffBasalName(self), self.kOffBasal );
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