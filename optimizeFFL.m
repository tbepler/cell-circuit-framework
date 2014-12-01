function optimizeFFL(plot)


set(0,'DefaultAxesFontSize', 20)
set(0,'DefaultTextFontSize', 20)
set(0,'DefaultLineLinewidth',3)

res = '-r900';
format = '-depsc';
ymin = 0;
ymax = 25;
%set(0,'DefaultFigureColor','w')

output = 'figures/ffl/optimize/edgeDetector_optimized';

if( ~exist( output, 'dir' ) )
    mkdir(output);
end

%diary( 'figures/ffl/optimize/edgeDetector_optimize_results.txt' );
%diary on

pulses = [ Pulse( 0, @(sys) sys ) ...
           Pulse( 20, @(sys) sys.ChangeInitialValue( 'X', 10 ) ) ...
           Pulse( 140, @(sys) sys.ChangeConstant( 'gamma_X', 10 ) ) ...
           Pulse( 200, @(sys) sys ) ];
       
    function o = obj( T )
        o = zeros( size(T) );
        for i = 1:numel(T)
            t = T(i);
            if t > 20 && t <= 30 
                o(i) = 10;
            end
        end
    end

lgd = { 'X', 'Y', 'Z', 'Y_{mRNA}', 'Z_{mRNA}'};

Ribo = 1000;

%y_unbindX z_unbindX z_unbindY x_coop y_coop
x0 = [ 5 5 10 2 4 ];
LB = [ 1 1 1 1 1 ];
UB = [ 50 50 50 8 8 ];

opts = psoptimset( 'Display', 'iter' );
[optParams, dist] = patternsearch( @distanceFromTarget, x0, ...
    [], [], [], [], LB, UB, [], opts );

P = optParams;
D = dist;

if(nargin > 0 && plot)
    
    [Y,T] = simulate( P );
    
    f = figure();
    hold on
    plot(T, Y);
    plot(T, obj(T), 'k');
    legend( [lgd, {'Objective'} ] );
    ylim([ymin ymax])
    xlabel('Time');
    ylabel('Concentration')
    hold off
    
    name = [ output ];
    for k = 1:length(params)
        name = [ name '_' num2str(params(k)) ];
    end
    
    print(f, name, res, format);
    
end

    function dist = distanceFromTarget( params )
        [Y,T] = simulate( params );
        dist = distance( Y(:,3), T );
    end

    function d = distance( Y, T )
        %d = sum( ( Y - obj(T) ) .^ 2 );
        objHigh = T > 20 & T <= 30;
        pulse = T > 20 & T <= 140;
        p = sum( Y(objHigh) ) / sum(objHigh);
        np = sum( Y(pulse & ~objHigh) ) / (sum(pulse & ~objHigh) );
        d = -(p / np);
    end

    function [Y,T] =  simulate( params )
        
        y_unbindX = params(1);
        z_unbindX = params(2);
        z_unbindY = params(3);
        x_coop = params(4);
        y_coop = params(5);
        
        [sys, x, y, z, y_mrna, z_mrna ] = createIncoherentFFL( ...
            100 ... %total RNAP
            , Ribo ... %total Ribo
            , 10 ... % k_txn
            , 1 ... % k_tln
            , 0 ... % x_deg
            , 0 ... % x_init
            , 2 ... % y_deg
            , 0 ... % y_init
            , 1 ... % z_deg
            , 0 ... % z_init
            , 5 ... % y_rna_deg
            , 1 ... % y_ribo_on
            , 50 ... % y_ribo_off
            , 0 ... % y_rna_init
            , 5 ... % z_rna_deg
            , 1 ... % z_ribo_on
            , 50 ... % z_ribo_off
            , 0 ... % z_rna_init
            , 1 ... % y_rnap_on
            , 10 ... % y_rnap_off
            , 1 ... y_bindX
            , y_unbindX ... y_unbindX
            , 1 ... y_copy_number
            , 1 ... z_rnap_on
            , 10 ... z_rnap_off
            , 1 ... z_bindX
            , z_unbindX ... z_unbindX
            , 1 ... z_bindY
            , z_unbindY ... z_unbindY
            , 1 ... zy_bindX
            , z_unbindX ... zy_unbindX
            , 1 ... zx_bindY
            , z_unbindY ... zx_unbindY
            , 1 ... z_copy_number
            , x_coop ... x_coop
            , y_coop ... y_coop
            );
        
        [T,Y] = sys.run_pulses( pulses );
        
        indices = [ sys.indexOf(x) sys.indexOf(y) sys.indexOf(z) ...
           sys.indexOf(y_mrna) sys.indexOf(z_mrna) ...
           ];
       
        Y = Y(:, indices);
        
    end

end