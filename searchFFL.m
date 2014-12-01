function searchFFL()

set(0,'DefaultAxesFontSize', 16)
set(0,'DefaultTextFontSize', 16)
set(0,'DefaultLineLinewidth',2)

res = '-r900';
format = '-depsc';
ymin = 0;
ymax = 25;
%set(0,'DefaultFigureColor','w')

output = 'figures/ffl/edgeDetector_gridSearch';

if( ~exist( output, 'dir' ) )
    mkdir(output);
end

pulses = [ Pulse( 0, 'X', 0 ) ...
           Pulse( 20, 'X', 10 ) ...
           Pulse( 140, 'X', 0 ) ...
           Pulse( 141, 'X', 0 ) ...
           Pulse( 142, 'X', 0 ) ...
           Pulse( 143, 'X', 0 ) ...
           Pulse( 144, 'X', 0 ) ...
           Pulse( 145, 'X', 0 )...
           Pulse( 146, 'X', 0 ) ...
           Pulse( 147, 'X', 0 ) ...
           Pulse( 148, 'X', 0 ) ...
           Pulse( 149, 'X', 0 ) ...
           Pulse( 150, 'X', 0 ) ...
           Pulse( 200, 'X', 0 ) ];
       
       
       
    function o = obj( T )
        o = zeros( size(T) );
        for i = 1:numel(T)
            t = T(i);
            if t > 20 && t <= 30 
                o(i) = 10;
            end
        end
    end

%Zt = [ zeros( 20, 1 ); repmat( 10, 10, 1 ); zeros( 170, 1 ) ];
       
lgd = { 'X', 'Y', 'Z', 'Y_{mRNA}', 'Z_{mRNA}'};

Ribo = 1000;

unbind = [ 5 10 20 ];
coop = [ 1 2 4 8 ];

[y_unbindXGrid, z_unbindXGrid, z_unbindYGrid, x_coopGrid, y_coopGrid ] = ndgrid( unbind, unbind, unbind, coop, coop );
V = gridSearch( @distanceFromTarget, y_unbindXGrid, z_unbindXGrid ...
    , z_unbindYGrid, x_coopGrid, y_coopGrid );

[m,I] = min(V);

y_unbindXGrid(I)
z_unbindXGrid(I)
z_unbindYGrid(I)
x_coopGrid(I)
y_coopGrid(I)
m

    function dist = distanceFromTarget( params )
        [Y,T] = simulate( params );
        dist = distance( Y(3), obj(T) );
    end

    function d = distance( X, Xt )
        d = sum( ( X - Xt ) .^ 2 );
    end

    function [Y,T] =  simulate( params )
        
        y_unbindX = params(1);
        z_unbindX = params(2);
        z_unbindY = params(3);
        x_coop = params(4);
        y_coop = params(5);
        msg = 'Trying params: ';
        for k = 1:length(params)
            msg = [msg ' ' num2str(params(k))];
        end
        disp( msg );
        
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
        
        %f = figure( 'Visible', 'off' );
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

end