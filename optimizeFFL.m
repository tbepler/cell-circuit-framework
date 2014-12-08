function [P,D] = optimizeFFL(plotResults)


set(0,'DefaultAxesFontSize', 20)
set(0,'DefaultTextFontSize', 20)
set(0,'DefaultLineLinewidth',2.5)

res = '-r900';
format = '-depsc';
ymin = 0;
ymax = 25;
%set(0,'DefaultFigureColor','w')

output = 'figures/ffl/optimize/';

if( ~exist( output, 'dir' ) )
    mkdir(output);
end

%diary( 'figures/ffl/optimize/edgeDetector_optimize_results.txt' );
%diary on

X_high_start = 20;
X_high_end = 140;
Z_duration_obj = 10;

pulses = [ Pulse( 0, @(sys) sys ) ...
           Pulse( X_high_start, @(sys) sys.ChangeInitialValue( 'X', 10 ) ) ...
           Pulse( X_high_end, @(sys) sys.ChangeConstant( 'gamma_X', 10 ) ) ...
           Pulse( 200, @(sys) sys ) ];
       
    function o = obj( T )
        o = zeros( size(T) );
        for i = 1:numel(T)
            t = T(i);
            if t > X_high_start && t <= X_high_start + Z_duration_obj
                o(i) = 15;
            end
        end
    end

lgd = { 'X', 'Y', 'Z', 'Y_{mRNA}', 'Z_{mRNA}'};

Ribo = 1000;

%y_deg z_deg y_ribo_off z_ribo_off z_unbindY y_coop
x0 = [ 2 2 50 50 10 4 ];
LB = [ 1 1 1 1 1 1 ];
UB = [ 100 100 100 100 100 5 ];

    function [C,Ceq] = nonlcon( x )
        C = 0;
        Ceq = ~( x(4) == floor( x(4) ) && x(5) == floor( x(5) ) );
    end

opts = psoptimset( 'Display', 'iter', 'UseParallel', true, 'TolFun', 1e-3 );
[optParams, dist] = patternsearch( @distanceFromTarget, x0, ...
    [], [], [], [], LB, UB, [], opts );

P = optParams;
D = dist;

if(nargin > 0 && plotResults)
    
    [Y,T] = simulate( P );
    
    f = figure();
    hold on
    plot(T,Y);
    plot(T, obj(T), 'k');
    legend( [lgd, {'Objective'} ] );
    ylim([ymin ymax])
    xlabel('Time');
    ylabel('Concentration')
    hold off
    
    name = [ output 'edgeDetector_optimized' ];
    for k = 1:length(P)
        name = [ name '_' num2str(P(k)) ];
    end
    name = [ name '.eps' ];
    
    print(f, name, res, format);
    
end

    function dist = distanceFromTarget( params )
        [Y,T] = simulate( params );
        dist = distance( Y(:,3), T );
    end

    function d = distance( Y, T )
        %d = sum( ( Y - obj(T) ) .^ 2 );
        objHigh = T > X_high_start & T <= X_high_start + Z_duration_obj;
        pulse = T > X_high_start & T <= X_high_end;
        p = sum( Y(objHigh) ) / sum(objHigh);
        np = sum( Y(pulse & ~objHigh) ) / (sum(pulse & ~objHigh) );
        bg = sum( Y(~pulse) ) / sum(~pulse) ;
        d = - p + np + abs( np - bg );
    end

    function [Y,T] =  simulate( params )
        
        y_deg = params(1);
        z_deg = params(2);
        y_ribo_off = params(3);
        z_ribo_off = params(4);
        z_unbindY = params(5);
        y_coop = round( params(6) ); %cooperitivity must be an integer
        
        [sys, x, y, z, y_mrna, z_mrna ] = createIncoherentFFL( ...
            100 ... %total RNAP
            , Ribo ... %total Ribo
            , 10 ... % k_txn
            , 1 ... % k_tln
            , 0 ... % x_deg
            , 0 ... % x_init
            , y_deg ... % y_deg
            , 0 ... % y_init
            , z_deg ... % z_deg
            , 0 ... % z_init
            , 5 ... % y_rna_deg
            , 1 ... % y_ribo_on
            , y_ribo_off ... % y_ribo_off
            , 0 ... % y_rna_init
            , 5 ... % z_rna_deg
            , 1 ... % z_ribo_on
            , z_ribo_off ... % z_ribo_off
            , 0 ... % z_rna_init
            , 1 ... % y_rnap_on
            , 10 ... % y_rnap_off
            , 1 ... y_bindX
            , 2 ... y_unbindX
            , 1 ... y_copy_number
            , 1 ... z_rnap_on
            , 10 ... z_rnap_off
            , 1 ... z_bindX
            , 2 ... z_unbindX
            , 1 ... z_bindY
            , z_unbindY ... z_unbindY
            , 1 ... zy_bindX
            , 2 ... zy_unbindX
            , 1 ... zx_bindY
            , z_unbindY ... zx_unbindY
            , 1 ... z_copy_number
            , 1 ... x_coop
            , y_coop ... y_coop
            );
        
        [T,Y] = sys.run_pulses( pulses );
        
        indices = [ sys.indexOf(x) sys.indexOf(y) sys.indexOf(z) ...
           sys.indexOf(y_mrna) sys.indexOf(z_mrna) ...
           ];
       
        Y = Y(:, indices);
        
    end

end