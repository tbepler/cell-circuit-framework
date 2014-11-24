function testSimpleResourceModel()

set(0,'DefaultAxesFontSize', 16)
set(0,'DefaultTextFontSize', 16)
set(0,'DefaultLineLinewidth',2)

res = '-r900';
format = '-depsc';
ymin = 0;
ymax = 40;
%set(0,'DefaultFigureColor','w')

output = 'figures/simple_model_test/';

if( ~exist( output, 'dir' ) )
    mkdir(output);
end

sys = SimpleResourceBioSystem( 8, 2, 2, 2 );
p = sys.addProtein( 'Protein', 0.1, 0 );
m = sys.addmRNA( 'mRNA', p, 0.05, 1, 0.5, 0 );
%addGene( name, mRNA, promoterStates, rnapOn, rnapOff, transitionMatrix, copy# )
g = sys.addGene( 'Gene', m, { [] }, [ 0.1 ], [ 1 ], zeros( 1, 1 ), 1 );

indices = [ sys.indexOf(p) sys.indexOf(m) sys.indexOf(g) ...
    sys.indexOfRNAP() sys.indexOfRibo() ];
lgd = { 'Protein', 'mRNA', 'Gene', 'RNAP', 'Ribosome' };

t = 0:100;
[T,Y] = sys.run( t );

%f = figure( 'Visible', 'off' );
f = figure();
hold on
plot(T, Y(:,indices));
legend( lgd );
ylim([ymin ymax]) 
xlabel('Time');
ylabel('Concentration')
hold off

name = [ output 'twoRibo_oneGene' ];
print(f, name, res, format);

%m.kOnRibo = 0.1;
%g.kOffBasal = 1.5;
%sys.setInitialValue( g, 5 );
sys.setRiboConcentration( 1 );


t = 0:100;
[T,Y] = sys.run( t );

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

sys.setInitialValue( g, 5 );

t = 0:100;
[T,Y] = sys.run( t );

%f = figure( 'Visible', 'off' );
f = figure();
hold on
plot(T, Y(:,indices));
legend( lgd );
ylim([ymin ymax]) 
xlabel('Time');
ylabel('Concentration')
hold off

name = [ output 'oneRibo_fiveGene' ];
print(f, name, res, format);