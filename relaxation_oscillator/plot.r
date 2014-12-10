library( ggplot2 )
library( reshape2 )

#define plotting function
plotFnc <- function( table ){

    fig <- ggplot( melt( table, "Time" ) , aes( x=Time,y=value,group=variable,colour=variable))+geom_line() + labs(y="Molecules",title="Simulation")
    return( fig )

}

args <- commandArgs( trailingOnly = TRUE )

print( "Plotting simulation..." )

if(length( args ) > 1){
    table <- read.table( args[1], header=TRUE )
    name <- args[2];
}else{
    table <- read.table( file( 'stdin' ), header=TRUE )
    name <- args[1];
}

#split file root and extension
split <- strsplit( name, "\\." )
root <- split[[1]][1]
ext <- split[[1]][2]

all_name <- paste( c( root, "_all.", ext ), collapse="" )
mrna_name <- paste( c( root, "_mrna.", ext ), collapse="" )

fig <- plotFnc( table ) 

ggsave(file=all_name,plot=fig,width=9,height=6,dpi=200)

keep <- c( "Time", "mRNA_X", "mRNA_Y" )
fig <- plotFnc( table[keep] )

ggsave(file=mrna_name,plot=fig,width=9,height=6,dpi=200 )


