function simulateRO(arg)
disp(['Starting simulation. Arg: ' arg ] );
WL = 200;
AMP = 10000;
timespan = 1 : 600;
plotevery = 12;

set(0,'DefaultAxesFontSize', 16)
set(0,'DefaultTextFontSize', 16)
set(0,'DefaultLineLinewidth',2)

res = '-r900';
format = '-depsc';
%set(0,'DefaultFigureColor','w')

output = 'figures/relaxation-oscillator/';

if( ~exist( output, 'dir' ) )
    mkdir(output);
end

switch arg
    case 'default'
        disp( 'Building default system...' );
        
%         params = Parameters.RELAXATION_OSCILLATOR_DEFAULT_PARAMETERS;
%         [ sys, pX, pY, mRNA_X, mRNA_Y, gene_X, gene_Y ] = createRelaxationOscillator( params{:} );
%         indices = [ sys.indexOf(pX) sys.indexOf(pY) sys.indexOf(mRNA_X) ...
%             sys.indexOf(mRNA_Y) ...
%             ];
%         lgd = { 'X', 'Y', 'X_{mRNA}', 'Y_{mRNA}'};
        disp('Building system');
        sys = RelaxationOscillator( RelaxationOscillator.DEFAULTS{:} );
        
        disp( 'Running simulation...' );
        
        %[T,Y] = sys.run( timespan );
        [ T, Y ] = sys.simulate( timespan );
        
        f = figure();
        hold on
        plot(T, Y(:,1:2));
        plot(T, objective(T), 'k');
        lgd = { 'X', 'Y','Objective'};
        legend( lgd );
        %ylim([ymin ymax])
        xlabel('Time');
        ylabel('Concentration')
        hold off
        
        name = [ output n '_defaults' ];
        print(f, name, res, format);
        
    case 'optim'
        disp( 'Optimizing parameters...');
        x0 = [
            Parameters.DEFAULT_PROT_DEG
            Parameters.MAX_PROT_DEG
            Parameters.DEFAULT_TF_COOP
            Parameters.DEFAULT_TF_COOP
            0.9
            0.2
            0.2
            Parameters.DEFAULT_RNA_DEG
            Parameters.DEFAULT_RIBO_INIT
            Parameters.DEFAULT_RNA_DEG
            Parameters.DEFAULT_RIBO_INIT
            ];
        LB = [
            Parameters.MIN_PROT_DEG
            Parameters.MIN_PROT_DEG
            Parameters.MIN_TF_COOP
            Parameters.MIN_TF_COOP
            0
            0
            0
            Parameters.MIN_RNA_DEG
            Parameters.MIN_RIBO_INIT
            Parameters.MIN_RNA_DEG
            Parameters.MIN_RIBO_INIT
        ];
        UB = [ 
            Parameters.MAX_PROT_DEG
            Parameters.MAX_PROT_DEG
            Parameters.MAX_TF_COOP
            Parameters.MAX_TF_COOP
            1
            1
            1
            Parameters.MAX_RNA_DEG
            Parameters.MAX_RIBO_INIT
            Parameters.MAX_RNA_DEG
            Parameters.MAX_RIBO_INIT
        ];
        
        evals = 0;
        opts = psoptimset( 'Display', 'iter', 'UseParallel', false, 'TolFun', 1e-3 );
        [optParams, dist] = patternsearch( @evaluate, x0, ...
            [], [], [], [], LB, UB, [], opts );
%         opts = optimoptions( 'fmincon', 'Display', 'iter', 'TolFun', 1e-3 );
%         optParams = fmincon( @evaluate, x0, ...
%             [], [], [], [], LB, UB, [], opts );
        params = makeParameters( optParams );
        [ sys, pX, pY, mRNA_X, mRNA_Y, gene_X, gene_Y ] = createRelaxationOscillator( params{:} );
        
        indices = [ sys.indexOf(pX) sys.indexOf(pY) sys.indexOf(mRNA_X) ...
            sys.indexOf(mRNA_Y) ...
            ];
        
        [T,Y] = sys.run( timespan );
        plotSys( T, Y(:,indices), 'optimal' );
        
        
end

    function p = makeParameters( vec )
        x_deg = vec(1);
        y_deg = vec(2);
        x_coop = round(vec(3));
        y_coop = round(vec(4));
        y_affY = vec(5);
        y_affX = vec(6);
        x_aff = vec(7);
        x_rna_deg = vec(8);
        x_ribo_init = vec(9);
        y_rna_deg = vec(10);
        y_ribo_init = vec(11);
        
        p = Parameters.RELAXATION_OSCILLATOR_DEFAULT_PARAMETERS;
        p{5} = x_deg;
        p{7} = y_deg;
        p{36} = x_coop;
        p{37} = y_coop;
        p{33} = Parameters.MIN_TF_ON + ( Parameters.MAX_TF_ON - Parameters.MIN_TF_ON ) * y_affY;
        p{34} = Parameters.MAX_TF_OFF - ( Parameters.MAX_TF_OFF - Parameters.MIN_TF_OFF ) * y_affY;
        p{31} = Parameters.MIN_TF_ON + ( Parameters.MAX_TF_ON - Parameters.MIN_TF_ON ) * y_affX;
        p{32} = Parameters.MAX_TF_OFF - ( Parameters.MAX_TF_OFF - Parameters.MIN_TF_OFF ) * y_affX;
        p{22} = Parameters.MIN_TF_ON + ( Parameters.MAX_TF_ON - Parameters.MIN_TF_ON ) * x_aff;
        p{23} = Parameters.MAX_TF_OFF - ( Parameters.MAX_TF_OFF - Parameters.MIN_TF_OFF ) * x_aff;
        p{9} = x_rna_deg;
        p{12} = x_ribo_init;
        p{14} = y_rna_deg;
        p{17} = y_ribo_init;
    end

    function plotSys( T, Y, n )
        %f = figure( 'Visible', 'off' );
        f = figure();
        hold on
        plot(T, Y(:,1:2));
        plot(T, objective(T), 'k');
        lgd = { 'X', 'Y','Objective'};
        legend( lgd );
        %ylim([ymin ymax])
        xlabel('Time');
        ylabel('Concentration')
        hold off
        
        name = [ output n '_proteins' ];
        print(f, name, res, format);
        
        f = figure();
        hold on
        plot(T, Y(:,3:4));
        %plot(T, objective(T), 'k');
        lgd = { 'X_{mRNA}', 'Y_{mRNA}' ...
            ...,'Objective' 
            };
        legend( lgd );
        %ylim([ymin ymax])
        xlabel('Time');
        ylabel('Concentration')
        hold off
        
        name = [ output n '_mrnas' ];
        print(f, name, res, format);
    end
    
    function F = objective( T )
        F = AMP .* cos( 2*pi.*T ./ WL ) ./2 + AMP/2;
    end

    function o = difference( X, Y )
        o = sum( abs( X - Y ) ) / length( X );
    end

    function score = evaluate( vec )
        evals = evals + 1;
        p = makeParameters( vec );
        
        [ sys, pX, pY, mRNA_X, mRNA_Y, gene_X, gene_Y ] = createRelaxationOscillator( p{:} );
        
        indices = [ sys.indexOf(pX) sys.indexOf(pY) sys.indexOf(mRNA_X) ...
            sys.indexOf(mRNA_Y) ...
            ];
        
        [T,Y] = sys.run( timespan );
        score = difference( objective(T), Y(:,sys.indexOf(pY) ) );
        
        d = evals / plotevery;
        if d == round( d )
            p'
            plotSys( T, Y(:,indices), [ 'evaluation' num2str(evals) ] );
        end
        
        
    end




% %m.kOnRibo = 0.1;
% %g.kOffBasal = 1.5;
% %sys.setInitialValue( g, 5 );
% sys.setRiboConcentration( 500 );
%
%
% %t = 0:100;
% [T,Y] = sys.run_pulses( pulses );
%
% %f = figure( 'Visible', 'off' );
% f = figure();
% hold on
% plot(T, Y(:,indices));
% legend( lgd );
% ylim([ymin ymax])
% xlabel('Time');
% ylabel('Concentration')
% hold off
%
% name = [ output 'oneRibo_oneGene' ];
% print(f, name, res, format);
%
% sys.setRiboConcentration( 250 );
%
% %t = 0:100;
% [T,Y] = sys.run_pulses( pulses );
%
% %f = figure( 'Visible', 'off' );
% f = figure();
% hold on
% plot(T, Y(:,indices));
% legend( lgd );
% ylim([ymin ymax])
% xlabel('Time');
% ylabel('Concentration')
% hold off
%
% name = [ output 'oneRibo_threeGene' ];
% print(f, name, res, format);

end