library( ggplot2 )
library( reshape2 )

args <- commandArgs( trailingOnly = TRUE )

print( "Hello world" )

if(length( args ) > 1){
    table <- read.table( args[1], header=TRUE )
    name <- args[2];
}else{
    table <- read.table( file( 'stdin' ), header=TRUE )
    name <- args[1];
}

fig <- ggplot( melt( table, "Time" ) , aes( x=Time,y=value,group=variable,colour=variable))+geom_line() +
labs(y="Molecules",title="Simulation")

ggsave(file=name,plot=fig,width=9,height=6,dpi=200)


