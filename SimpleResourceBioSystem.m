classdef SimpleResourceBioSystem < BioSystem
    
    properties(Constant = true)
        RNAP_NAME = 'RNAP';
        RIBO_NAME = 'Ribo';
        K_TXN_NAME = 'k_txn';
        K_TLN_NAME = 'k_tln';
    end
    
    properties(SetAccess = immutable)
        RNAP
        Ribo
    end
    
    properties(Access = private)
        proteins
        mRNAs
        genes
    end
    
    methods
        function obj = SimpleResourceBioSystem( concRNAP, concRibo ...
                , kTxn, kTln )
            obj = obj@BioSystem();
            obj.RNAP = Compositor( obj.RNAP_NAME, concRNAP );
            obj.Ribo = Compositor( obj.RIBO_NAME, concRibo );
            
            obj.AddCompositor( obj.RNAP );
            obj.AddCompositor( obj.Ribo );
            
            obj.AddConstant( obj.K_TXN_NAME, kTxn ); %CURRENTLY IGNORED
            obj.AddConstant( obj.K_TLN_NAME, kTln ); %CURRENTLY IGNORED
            
            obj.proteins = containers.Map;
            obj.mRNAs = containers.Map;
            obj.genes = containers.Map;
            
        end
        
        function self = setTxnRate( self, kTxn )
            ChangeConstant( self, self.K_TXN_NAME, kTxn );
        end
        
        function self = setTlnRate( self, kTln )
            ChangeConstant( self, self.K_TLN_NAME, kTln );
        end
        
        function self = setRNAPConcentration( self, conc )
            self.ChangeInitialValue( self.RNAP_NAME, conc );
        end
        
        function self = setRiboConcentration( self, conc )
            self.ChangeInitialValue( self.RIBO_NAME, conc );
        end
        
        function self = setInitialValue( self, obj, val )
            self.ChangeInitialValue( obj.name, val );
        end
        
        function comp = getCompositor( self, obj )
            i = self.map_compositors( obj.name );
            comp = self.compositors(i);
        end
        
        function i = indexOf( self, obj )
            i = self.map_compositors( obj.name );
        end;
        
        function i = indexOfRibo( self )
            i = self.map_compositors( self.RIBO_NAME );
        end
        
        function i = indexOfRNAP( self )
            i = self.map_compositors( self.RNAP_NAME );
        end
        
        function p = addProtein( self, varargin )
            assert( length(varargin) == 2 || length(varargin) == 3 ...
                , 'Requires %d or %d args', 2, 3 );
            if( length(varargin) == 2 )
                p = varargin{1};
                conc = varargin{2};
            elseif( length(varargin) == 3 )
                p = Protein( varargin{1}, varargin{2} );
                conc = varargin{3};
            end
            assert( ~isKey( self.proteins, p.name ) ...
                    , 'System already contains protein named "%s"' ...
                    , p.name);
            self.proteins(p.name) = p;
            
            p.accept( self, conc );
            
            addlistener( p, 'kDeg', 'PostSet', @self.update );
        end     
        
        function m = addmRNA( self, varargin )
            assert( length(varargin) == 2 || length(varargin) == 7 ...
                , 'Requires %d or %d args', 2, 7 );
            if( length(varargin) == 2 )
                m = varargin{1};
                conc = varargin{2};
            elseif( length(varargin) == 7 )
                m = mRNA( varargin{1}, varargin{2}, varargin{3} ...
                    , varargin{4}, varargin{5}, varargin{6} );
                conc = varargin{7};
            end
            assert( ~isKey( self.mRNAs, m.name ) ...
                    , 'System already contains mRNA named "%s"' ...
                    , m.name);
            self.mRNAs(m.name) = m;
            
            m.accept( self, conc );
            
            addlistener( m, 'kDeg', 'PostSet', @self.update );
            addlistener( m, 'kOnRibo', 'PostSet', @self.update );
            addlistener( m, 'kOffRibo', 'PostSet', @self.update );
        end   
        
        function g = addGene( self, varargin )
            assert( length(varargin) == 2 || length(varargin) == 8 ...
                , 'Requires %d or %d args', 2, 8 );
            if( length(varargin) == 2 )
                g = varargin{1};
                conc = varargin{2};
            elseif( length(varargin) == 8 )
                %Gene( name, mRNA, promoterStates, rnapOn, rnapOff,
                %rnapInit,
                %transitionMatrix )
                g = Gene( varargin{1}, varargin{2}, varargin{3} ...
                    , varargin{4}, varargin{5}, varargin{6}, varargin{7} );
                conc = varargin{8};
            end
            assert( ~isKey( self.genes, g.name ) ...
                    , 'System already contains gene named "%s"' ...
                    , g.name);
            self.genes(g.name) = g;
            
            g.accept( self, conc );
            
            %addlistener( g, 'kOnBasal', 'PostSet', @self.update );
            %addlistener( g, 'kOffBasal', 'PostSet', @self.update );
        end   
        
        function update( self, ~, e )
            e.AffectedObject.update( self );
        end
        
        function delete( self )
            self.proteins = [];
            self.mRNAs = [];
            self.genes = [];
        end
        
    end
end


