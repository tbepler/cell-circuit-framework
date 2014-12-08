classdef RelaxationOscillator
    %RELAXATIONOSCILLATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( Constant )
        
        DEFAULT_K_D = Parameters.DEFAULT_TF_OFF / Parameters.DEFAULT_TF_ON;
        DEFAULT_K_MRNAP = (Parameters.DEFAULT_RNAP_INIT + Parameters.DEFAULT_RNAP_OFF) ...
            / Parameters.DEFAULT_RNAP_ON ;
        DEFAULT_OCL = Parameters.DEFAULT_RNAP_INIT;
        
        
        DEFAULTS = {
            Parameters.DEFAULT_RNAP
            Parameters.DEFAULT_RIBO
            1
            1
            Parameters.DEFAULT_TF_COOP
            Parameters.DEFAULT_TF_COOP
            RelaxationOscillator.DEFAULT_K_D;
            RelaxationOscillator.DEFAULT_K_MRNAP;
            RelaxationOscillator.DEFAULT_OCL;
            Parameters.DEFAULT_RNA_DEG;
            RelaxationOscillator.DEFAULT_K_D;
            RelaxationOscillator.DEFAULT_K_D;
            RelaxationOscillator.DEFAULT_K_MRNAP;
            RelaxationOscillator.DEFAULT_OCL;
            RelaxationOscillator.DEFAULT_K_MRNAP;
            RelaxationOscillator.DEFAULT_OCL;
            Parameters.DEFAULT_RNA_DEG;
            Parameters.DEFAULT_RIBO_ON;
            Parameters.DEFAULT_RIBO_OFF;
            Parameters.DEFAULT_RIBO_INIT;
            Parameters.DEFAULT_RIBO_INIT;
            Parameters.DEFAULT_PROT_DEG;
            Parameters.DEFAULT_PROT_DEG;
        };
        
    end
    
    properties
        
        rateY;
        rateX;
        
    end
    
    methods
        
        function obj = RelaxationOscillator( RNAP, Ribo, G_X, G_Y, X_COOP ...
                , Y_COOP , K_DXY, K_MXRNAP, K_XOCL, K_RNAX_DEG, K_DYY, K_DYX ...
                , K_MYRNAP_BASAL, K_YOCL_BASAL, K_MYRNAP_ACTIVE, K_YOCL_ACTIVE ...
                , K_RNAY_DEG, K_ON_RIBO, K_OFF_RIBO, K_TLNX, K_TLNY, K_DEG_X, K_DEG_Y)
            
            syms X Y;
            disp( 'create mrna x ' );
            mRNA_X = mRNA_X_Eq( G_X, Y^Y_COOP, RNAP, K_DXY^Y_COOP, K_MXRNAP ... 
                , K_XOCL, K_RNAX_DEG );
            
            disp( 'create mrna y ' );
            mRNA_Y = mRNA_Y_Eq( G_Y, Y^Y_COOP, X^X_COOP, RNAP, K_DYY^Y_COOP ...
                , K_DYX^X_COOP, K_MYRNAP_BASAL ...
                , K_YOCL_BASAL, K_MYRNAP_ACTIVE, K_YOCL_ACTIVE, K_RNAY_DEG );
            
            
            disp( 'create rate x ' );
            symRateX = proteinRate( X, mRNA_X, Ribo, K_ON_RIBO, K_OFF_RIBO ...
                , K_TLNX, K_DEG_X );
            
            
            disp( 'create mrna y ' );
            symRateY = proteinRate( Y, mRNA_Y, Ribo, K_ON_RIBO, K_OFF_RIBO ...
                , K_TLNY, K_DEG_Y );
            
            
            disp( 'compile x ' );
            obj.rateX = matlabFunction( simplify( symRateX ) );
            obj.rateX
            
            disp( 'compile y' );
            obj.rateY = matlabFunction( simplify( symRateY ) );
            obj.rateY
            
        end
        
        function [T,Y] = simulate( self, tspan, initX, initY )
            y0 = [ initX ; initY ];
            [ T, Y ] = ode23s(@self.eval, tspan, y0);
        end
        
        function dy = eval( self, ~, y )
            dy = [
                self.rateX( y(1), y(2) );
                self.rateY( y(1), y(2) );
            ];
        end
        
    end
    
end

function r = proteinRate( P, mRNA, Ribo, k_on, k_off, k_cat, k_deg )

K_M = ( k_off + k_cat ) / k_on;
tlnComplex = Reactions.complexAtEq( mRNA, Ribo, K_M );
tlnRate = k_cat * tlnComplex;
degRate = k_deg * P;

r = tlnRate - degRate;

end

function rna = mRNA_X_Eq( G_X, Y, RNAP, K_DY, K_Mrnap, k_ocl, k_deg )

activeComplex = Reactions.complexAtEq( G_X, Y, K_DY );
rnaProd = Reactions.mrnaProduction( activeComplex, RNAP, K_Mrnap, k_ocl );
rna = Reactions.mrnaAtEq( [ rnaProd ], k_deg );

end

function rna = mRNA_Y_Eq( G_Y, Y, X, RNAP, K_DY, K_DX ...
    , K_Mbasal, basal_ocl, K_Mactive, active_ocl ...
    , k_deg )

boundY = Reactions.complexAtEq( G_Y, Y, K_DY );
boundX = Reactions.complexAtEq( G_Y, X, K_DX );
boundXY = Reactions.complexAtEq( boundY, X, K_DX );
naked = G_Y - boundY - boundX + boundXY;
boundYActive = boundY - boundXY;

boundYRNAProd = Reactions.mrnaProduction( boundYActive, RNAP, K_Mactive, active_ocl );

basalRNAProd = Reactions.mrnaProduction( naked, RNAP, K_Mbasal, basal_ocl );

rna = Reactions.mrnaAtEq( [ boundYRNAProd, basalRNAProd ], k_deg );

end