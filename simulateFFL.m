function simulateFFL()

set(0,'DefaultAxesFontSize', 16)
set(0,'DefaultTextFontSize', 16)
set(0,'DefaultLineLinewidth',2)

res = '-r900';
format = '-depsc';
ymin = 0;
ymax = 25;
%set(0,'DefaultFigureColor','w')

output = 'figures/ffl/edgeDetector_';

if( ~exist( output, 'dir' ) )
    mkdir(output);
end

disp( 'Building system...' );

sys = SimpleResourceBioSystem( 100, 1000, 10, 1 );
%addProtein( name, degredation rate, initial value )
x = sys.addProtein( 'X', 0, 0 );
y = sys.addProtein( 'Y', 2, 0);
z = sys.addProtein( 'Z', 1, 0);
%addmRNA( name, protein, degredation rate, ribo on rate, ribo off rate, initial value )
y_mrna = sys.addmRNA( 'mRNA_Y', y, 5, 1, 50, 0 );
z_mrna = sys.addmRNA( 'mRNA_Z', z, 5, 1, 50, 0 );
%addGene( name, mRNA, promoterStates, rnapOn, rnapOff, transitionMatrix, copy# )
y_gene = sys.addGene( 'Gene_Y', y_mrna, { [], [x,x] }, [ 0, 1 ], [ 0, 10 ], [ 0 , 1 ; 25, 0 ] , 1 );
z_states = { [], [x,x], [y, y, y, y], [x, x, y, y, y, y] };
%edge detector ffl
z_on = [ 0, 1, 0, 0 ];
z_off = [ 0, 10, 0, 0 ];
z_trans = [ 0 , 1, 1, 0 ; ...
            25 , 0, 0, 1 ; ...
            4^4, 0, 0, 1 ; ...
            0, 4^4, 25, 0 ];
z_gene = sys.addGene( 'Gene_Z', z_mrna, z_states, z_on, z_off, z_trans, 1 ); 

indices = [ sys.indexOf(x) sys.indexOf(y) sys.indexOf(z) ...
    sys.indexOf(y_mrna) sys.indexOf(z_mrna) ...
    ];
lgd = { 'X', 'Y', 'Z', 'Y_{mRNA}', 'Z_{mRNA}'};

disp( 'Running simulation...' );

%t = 0:100;
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
           %Pulse( 80, 'X', 0 ) ...
           %Pulse( 100, 'X', 0 ) ]; 
[T,Y] = sys.run_pulses( pulses );

%f = figure( 'Visible', 'off' );
f = figure();
hold on
plot(T, Y(:,indices));
legend( lgd );
ylim([ymin ymax]) 
xlabel('Time');
ylabel('Concentration')
hold off

name = [ output 'eightRibo_oneGene' ];
print(f, name, res, format);

%m.kOnRibo = 0.1;
%g.kOffBasal = 1.5;
%sys.setInitialValue( g, 5 );
sys.setRiboConcentration( 500 );


%t = 0:100;
[T,Y] = sys.run_pulses( pulses );

%f = figure( 'Visible', 'off' );
f = figure();
hold on
plot(T, Y(:,indices));
legend( lgd );
ylim([ymin ymax]) 
xlabel('Time');
ylabel('Concentration')
hold off

name = [ output 'oneRibo_oneGene' ];
print(f, name, res, format);

sys.setRiboConcentration( 250 );

%t = 0:100;
[T,Y] = sys.run_pulses( pulses );

%f = figure( 'Visible', 'off' );
f = figure();
hold on
plot(T, Y(:,indices));
legend( lgd );
ylim([ymin ymax]) 
xlabel('Time');
ylabel('Concentration')
hold off

name = [ output 'oneRibo_threeGene' ];
print(f, name, res, format);