function [sys, RNAP, Ribo] = setupSimpleResourceModel( concRNAP, concRibo )

sys = BioSystem();

RNAP = Compositor('RNAP', concRNAP);
Ribo = Compositor('Ribo', concRibo);

sys.AddCompositor( RNAP );
sys.AddCompositor( Ribo );
