#include<stdio.h>
#include<math.h>
#include<string.h>

#define RNAP_ECOLI  4600
#define RIBO_ECOLI  55000
#define TXN  1050
#define TLN  1050
#define MAX_DEG_PROTEIN  1.37e-2
#define MIN_DEG_PROTEIN  5e-4
#define MAX_DEG_RNA  2.31e-1
#define MIN_DEG_RNA  7.7e-3
#define RIBO_ON  1.5398
#define RIBO_OFF  0.9
#define MAX_RIBO_INIT  3.6
#define MIN_RIBO_INIT  1.2
#define MAX_RNAP_KM  1.2118e5
#define MIN_RNAP_KM  116.0585
#define MAX_RNAP_INIT  3.6
#define MIN_RNAP_INIT  1.2
#define MAX_TF_KD  6.6243e3
#define MIN_TF_KD  0.6624
#define MIN_TF_COOP  1
#define MAX_TF_COOP  10


typedef struct Constants{
    
    /* Protein constant */
    double deg_y;
    double deg_x;
    double coop_y;
    double coop_x;

    /* mRNA constants */
    double deg_mrna_y;
    double deg_mrna_x;

    double ribo_on_mrna_y;
    double ribo_off_mrna_y;
    double tln_mrna_y;

    double ribo_on_mrna_x;
    double ribo_off_mrna_x;
    double tln_mrna_x;

    /* Gene constants */
    double gene_x_kd_y;
    double gene_x_km_rnap;
    double gene_x_txn;

    double gene_y_kd_y;
    double gene_y_kd_x;
    double gene_y_basal_km_rnap;
    double gene_y_basal_txn;
    double gene_y_active_km_rnap;
    double gene_y_active_txn;

} Constants;

typedef struct Species{
    
    /* Resources */
    double rnap;
    double ribo;

    /* Proteins */
    double y;
    double x;
    
    /* mRNA */
    double mrna_y;
    double mrna_x;

    /* genes */
    double gene_y;
    double gene_x;
    
} Species;


static inline double complex( double a, double b, double kd ){
/*

Solution to:
A * B / ( Kd + B ) == C
where
A + C == A_total
B + C == B_total

C == ( A + B + Kd  - sqrt( A^2 - 2*A*B + 2*A*Kd + B^2 + 2*B*Kd + Kd^2 ) ) / 2

*/

return ( a + b + kd - sqrt( pow( a, 2 ) - 2*a*b + 2*a*kd + pow( b, 2 ) + 2*b*kd + pow( kd, 2 ) ) ) / 2.0; 

}

static inline double proteinRate( double prot, double mrna, double ribo, double deg, double tln, double ribo_on, double
ribo_off ){
    
    // dP/dt = tln * [mRNA-Ribo] - deg * P
    double km = ( ribo_off + tln ) / ribo_on;
#ifdef DEBUG
    fprintf( stderr, "ribo=%lf, mrna=%lf, km=%lf\n", ribo, mrna, km );
    fprintf( stderr, " %lf * %lf - %lf * %lf = %lf\n", tln, complex( mrna, ribo, km ), deg, prot,
    tln*complex(mrna,ribo,km) - deg*prot );
#endif
    return tln * complex( mrna, ribo, km ) - deg * prot;

}

/*  */
static inline void proteinDynamics( const Constants* c, const Species* cur, Species* next, double dt ){
    
    double dy = proteinRate( cur->y, cur->mrna_y, cur->ribo, c->deg_y, c->tln_mrna_y, c->ribo_on_mrna_y, c->ribo_off_mrna_y ) *
    dt;

    double dx = proteinRate( cur->x, cur->mrna_x, cur->ribo, c->deg_x, c->tln_mrna_x, c->ribo_on_mrna_x, c->ribo_off_mrna_x ) *
    dt;
    
    next->y = next->y + dy;
    next->x = next->x + dx;

}

static inline double mRNA_X_Production( double gene, double rnap, double tf, double coop, double kd, double km, double txn ){

    // dR/dt = txn * [ gene-tf-rnap]
    double gene_tf = complex( gene, pow( tf, coop ), pow( kd, coop ) );
    double gene_tf_rnap = complex( gene_tf, rnap, km );
 /*   fprintf( stderr, "tf=%lf, gene_tf=%lf, gene_tf_rnap=%lf, txn*gene_tf_rnap=%lf\n",tf, gene_tf, gene_tf_rnap,
    txn*gene_tf_rnap ); */
    return txn * gene_tf_rnap;

}

static inline double mRNA_Y_Production( double gene, double rnap, double activator, double activator_coop, double
activator_kd, double repressor, double repressor_coop, double repressor_kd, double km_basal, double txn_basal, double
km_active, double txn_active ){

    // dR/dt = txn_basal * [naked promoter] + txn_active*[promoter-y]
    double gene_activator = complex( gene, pow( activator, activator_coop ), pow( activator_kd, activator_coop ) );
    double gene_repressor = complex( gene, pow( repressor, repressor_coop ), pow( repressor_kd, repressor_coop ) );
    double gene_both = complex( gene_activator, pow( repressor, repressor_coop ), pow( repressor_kd, repressor_coop ) );
    double promoter_naked = gene - gene_activator - gene_repressor + gene_both;
    double promoter_active = gene_activator - gene_both;
    
    double promoter_naked_rnap = complex( promoter_naked, rnap, km_basal );
    double promoter_active_rnap = complex( promoter_active, rnap, km_active );

    return txn_basal * promoter_naked_rnap + txn_active * promoter_active_rnap;

}

static inline void mRNADynamics( const Constants* c, const Species* cur, Species* next, double dt ){

    double dmrna_x = ( mRNA_X_Production( cur->gene_x, cur->rnap, cur->y, c->coop_y, c->gene_x_kd_y, c->gene_x_km_rnap,
    c->gene_x_txn ) - c->deg_mrna_x * cur->mrna_x ) * dt;

    double dmrna_y = ( mRNA_Y_Production( cur->gene_y, cur->rnap, cur->y, c->coop_y, c->gene_y_kd_y, cur->x, c->coop_x,
    c->gene_y_kd_x, c->gene_y_basal_km_rnap, c->gene_y_basal_txn, c->gene_y_active_km_rnap, c->gene_y_active_txn ) -
    c->deg_mrna_y * cur->mrna_y ) * dt;

    next->mrna_x += dmrna_x;
    next->mrna_y += dmrna_y;

}

static inline void printHeader( FILE* out ){
    fprintf( out, "Time X Y mRNA_X mRNA_Y\n" );
}

static inline void report( FILE* out, const Species* cur, double time ){
    fprintf( out, "%lf %lf %lf %lf %lf\n", time, cur->x, cur->y, cur->mrna_x, cur->mrna_y );
}

void simulate( const Constants* c, const Species* init, double tend, double tstep, FILE* out, int reportIters ){
    
    printHeader( out );
    double time = 0;
    unsigned long iters = 0;
    Species cur;
    memcpy( &cur, init, sizeof( cur ) );
    Species next;
    memcpy( &next, init, sizeof( next ) );
    report( out, &cur, time );
    while( time < tend ){
        time += tstep;
        mRNADynamics( c, &cur, &next, tstep );
        proteinDynamics( c, &cur, &next, tstep );
        memcpy( &cur, &next, sizeof( cur ) );
        if( ++iters % reportIters == 0 ){
            report( out, &cur, time );
        }
    }
}

static const Constants initParams = {
    //1.67e-3, // deg_y
    MAX_DEG_PROTEIN,
    MAX_DEG_PROTEIN, // deg_x
    2, //coop_y
    4, //coop_x
    3.47e-2, //deg_mrna_y
    3.47e-2, //deg_mrna_x
    RIBO_ON, //ribo_on_mrna_y
    RIBO_OFF, //ribo_off_mrna_y
    (2.4 * TLN)/(2.4+TLN), //tln_mrna_y
    RIBO_ON, //ribo_on_mrna_x
    RIBO_OFF, //ribo_off_mrna_x
    (2.4 * TLN)/(2.4+TLN), //tln_mrna_x
    3500, //gene_x_kd_y
    2000, //gene_x_km_rnap
    (2.8 * TXN)/(2.8+TXN), //gene_x_txn
    1000, //gene_y_kd_y
    4000, //gene_y_kd_x
    3000, //gene_y_basal_km_rnap
    (1.2 * TXN)/(1.2+TXN), //gene_y_basal_txn
    1000, //gene_y_active_km_rnap
    (3.0 * TXN)/(3.0+TXN) // gene_y_active_txn
};

static const Species initValues = {
    RNAP_ECOLI, //rnap
    RIBO_ECOLI, //ribo
    0, //y
    0, //x
    0, //mrna_y
    0, //mrna_x
    1, //gene_y
    1 //gene_x
};

typedef enum Mode { SIMULATE, OPTIMIZE } Mode;

void parse( FILE* file, Constants* consts, Constants* minVals, Constants* maxVals, Species* initVals ){

    char tag[10];
    double val;
    double min;
    double max;

    while( !feof( file ) ) {
        int items = fscanf( file, "%s %lf %lf %lf", &tag, &val, &min, &max );
        if( items == 4 ){
            if( tag == "DEG_Y" ){
                consts->deg_y = val;
                minVals->deg_y = min;
                maxVals->deg_y = max;
            }
        }

    }

}

int main( int argc, char* argv[] ){
    
    FILE* out = stdout;
    if( argc > 1 ){
        out = fopen( argv[1], "r" );
        if( out == NULL ){
            fprintf( stderr, "Unable to open output file: %s\n", argv[1] );
            return 1;
        }
    }

    static const double TEND = 1000;
    static const double TSTEP = 0.01;
    static const int REPORT_ITERS = 100;

    fprintf( stderr, "Running simulation with TEND = %lf, TSTEP = %lf, REPORT_ITERS = %d\n", TEND, TSTEP, REPORT_ITERS
    );

    simulate( &initParams, &initValues, TEND, TSTEP, out, REPORT_ITERS );


}

