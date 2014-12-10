#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<getopt.h>
#include<assert.h>
#include"c_utils/file.h"
#include"c_utils/array.h"
#include"c_utils/optimize.h"

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

#define DEFAULT_TEND 2000
#define DEFAULT_TSTEP 0.1
#define DEFAULT_REPORT_ITERS 10

static const char* DELIMS = " \t";

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


//only valid for 1:1 binding
static inline double complex( double a, double b, double kd ){
    /*

       Solution to:
       A * B / ( Kd + B ) == C
       where
       A + C == A_total
       B + C == B_total

       C == ( A + B + Kd  - sqrt( A^2 - 2*A*B + 2*A*Kd + B^2 + 2*B*Kd + Kd^2 ) ) / 2

     */

    double b_2 = pow( b, 2 );
    double a_2 = pow( a, 2 );
    double kd_2 = pow( kd, 2 );
    double val = ( a + b + kd - sqrt( a_2 - 2*a*b + 2*a*kd + b_2 + 2*b*kd + kd_2 ) ) / 2;
    
    if( val > a || val > b ){ //numeric instability has occurred
        return a < b ? a : b;
    }
    return val;
}

static inline double degredation( double x, double k_deg, double dt ){
    
    //integrate dx/dt = - k_deg * x over timestep dt
    //double xt = x * exp( -k_deg * dt );
    return x * ( 1 - exp( -k_deg * dt ) );
    //return k_deg * x * dt;

}

static inline double proteinRate( double prot, double mrna, double ribo, double tln, double ribo_on, double
        ribo_off ){

    // dP/dt = tln * [mRNA-Ribo] - deg * P
    double km = ( ribo_off + tln ) / ribo_on;
    return tln * complex( mrna, ribo, km );

}

/*  */
static inline void proteinDynamics( const Constants* c, const Species* cur, Species* next, double dt ){

    double prodY = proteinRate( cur->y, cur->mrna_y, cur->ribo, c->tln_mrna_y, c->ribo_on_mrna_y, c->ribo_off_mrna_y ) * dt;
    double degY = degredation( cur->y, c->deg_y, dt );
    double dy = prodY - degY;

    double prodX = proteinRate( cur->x, cur->mrna_x, cur->ribo, c->tln_mrna_x, c->ribo_on_mrna_x, c->ribo_off_mrna_x ) * dt;
    double degX = degredation( cur->x, c->deg_x, dt );
    double dx = prodX - degX;

    next->y = next->y + dy;
    next->x = next->x + dx;

 }

static inline double mRNA_X_Production( double gene, double rnap, double tf, double coop, double kd, double km, double txn ){

    //hill eqn. for binding of both tf and RNAP
    double tf_coop = pow( tf, coop );
    double kd_coop = pow( kd, coop );
    double gene_tf_rnap = gene * ( tf_coop / ( kd_coop + tf_coop ) ) *
        ( rnap / ( km + rnap ) );

   return txn * gene_tf_rnap;

}

static inline double mRNA_Y_Production( double gene, double rnap, double activator, double activator_coop, double
        activator_kd, double repressor, double repressor_coop, double repressor_kd, double km_basal, double txn_basal, double
        km_active, double txn_active ){

    //Hill eqns to describe promoter states
    double act_coop = pow( activator, activator_coop );
    double act_kd_coop = pow( activator_kd, activator_coop );
    double rps_coop = pow( repressor, repressor_coop );
    double rps_kd_coop = pow( repressor_kd, repressor_coop );

    double promoter_naked_rnap = gene * ( act_kd_coop / ( act_kd_coop + act_coop ) ) *
        ( rps_kd_coop / ( rps_kd_coop + rps_coop ) ) * ( rnap / ( km_basal + rnap ) );
    double promoter_active_rnap = gene * ( act_coop / ( act_kd_coop + act_coop ) ) * 
        ( rps_kd_coop / ( rps_kd_coop + rps_coop ) ) * ( rnap / ( km_active + rnap ) );

    return txn_basal * promoter_naked_rnap + txn_active * promoter_active_rnap;

}

static inline void mRNADynamics( const Constants* c, const Species* cur, Species* next, double dt ){

    double prodX = mRNA_X_Production( cur->gene_x, cur->rnap, cur->y, c->coop_y, c->gene_x_kd_y, c->gene_x_km_rnap, c->gene_x_txn ) * dt;
    double degX = degredation( cur->mrna_x, c->deg_mrna_x, dt);
    double dmrna_x = prodX - degX;

    double prodY = mRNA_Y_Production( cur->gene_y, cur->rnap, cur->y, c->coop_y, c->gene_y_kd_y, cur->x, c->coop_x,c->gene_y_kd_x, c->gene_y_basal_km_rnap, c->gene_y_basal_txn, c->gene_y_active_km_rnap, c->gene_y_active_txn ) * dt;
    double degY = degredation( cur->mrna_y, c->deg_mrna_y, dt );
    double dmrna_y = prodY - degY;

    next->mrna_x += dmrna_x;
    next->mrna_y += dmrna_y;
    //ensure values never less than zero
    /*if( next->mrna_x < 0 ){
        next->mrna_x = 0;
    }
    if( next->mrna_y < 0 ){
        next->mrna_y = 0;
    }*/

    if( cur->mrna_y > 0 && next->mrna_y < 0 ){
        fprintf( stderr, "mrna_y_n=%lf, mrna_y_{n+1}=%lf, prod=%lf, deg=%lf\n", cur->mrna_y, next->mrna_y, prodY, degY );
    }
    if( cur->mrna_x > 0 && next->mrna_x < 0 ){
        fprintf( stderr, "mrna_x_n=%lf, mrna_x_{n+1}=%lf, prod=%lf, deg=%lf\n", cur->mrna_x, next->mrna_x, prodX, degX );
    }

}

static inline void printHeader( FILE* out ){
    fprintf( out, "Time X Y mRNA_X mRNA_Y\n" );
}

static inline void report( FILE* out, const Species* cur, double time ){
    fprintf( out, "%lf %lf %lf %lf %lf\n", time, cur->x, cur->y, cur->mrna_x, cur->mrna_y );
}

typedef struct Simulation{

    Species* vals;
    double* time;
    size_t n;

} Simulation;

Simulation simulate( const Constants* c, const Species* init, double tend, double tstep, FILE* out, int reportIters ){

    if( out != NULL ){
        printHeader( out );
    }
    double time = 0;
    unsigned long iters = 0;
    Species cur;
    memcpy( &cur, init, sizeof( cur ) );
    Species next;
    memcpy( &next, init, sizeof( next ) );

    //initialize Simulation
    Simulation sim;
    size_t steps = tend / tstep + 1;
    size_t n = steps / reportIters + 1;
    sim.vals = malloc( n * sizeof( Species ) );
    sim.time = malloc( n * sizeof( double ) );
    sim.n = n;
    
    memcpy( sim.vals, &cur, sizeof( *sim.vals ) );
    sim.time[0] = time;
    if( out != NULL ){
        report( out, &cur, time );
    }
    
    size_t idx = 1;

    while( steps-- > 0 ){
        time += tstep;
        mRNADynamics( c, &cur, &next, tstep );
        proteinDynamics( c, &cur, &next, tstep );
        memcpy( &cur, &next, sizeof( cur ) );
        if( ++iters % reportIters == 0 ){
            assert( idx < n );
            memcpy( sim.vals + idx, &cur, sizeof( *(sim.vals + idx ) ) );
            sim.time[idx] = time;
            ++idx;
            if( out != NULL ){
                report( out, &cur, time );
            }

        }
    }

    return sim;

}

static Species initValues = {
    RNAP_ECOLI, //rnap
    RIBO_ECOLI, //ribo
    0, //y
    0, //x
    0, //mrna_y
    0, //mrna_x
    1, //gene_y
    1 //gene_x
};

inline void arrayToConstants( const double* array, Constants* consts ){

    consts->deg_y = array[0];
    consts->deg_x = array[1];
    consts->coop_y = array[2];
    consts->coop_x = array[3];
    consts->deg_mrna_y = array[4];
    consts->deg_mrna_x = array[5];
    consts->ribo_on_mrna_y = array[6];
    consts->ribo_off_mrna_y = array[7];
    consts->tln_mrna_y = array[8];
    consts->ribo_on_mrna_x = array[9];
    consts->ribo_off_mrna_x = array[10];
    consts->tln_mrna_x = array[11];
    consts->gene_x_kd_y = array[12];
    consts->gene_x_km_rnap = array[13];
    consts->gene_x_txn = array[14];
    consts->gene_y_kd_y = array[15];
    consts->gene_y_kd_x = array[16];
    consts->gene_y_basal_km_rnap = array[17];
    consts->gene_y_basal_txn = array[18];
    consts->gene_y_active_km_rnap = array[19];
    consts->gene_y_active_txn = array[20];

}

inline void constantsToArray( const Constants* consts, double* array ){

    array[0] = consts->deg_y;
    array[1] = consts->deg_x;
    array[2] = consts->coop_y;
    array[3] = consts->coop_x;
    array[4] = consts->deg_mrna_y;
    array[5] = consts->deg_mrna_x;
    array[6] = consts->ribo_on_mrna_y;
    array[7] = consts->ribo_off_mrna_y;
    array[8] = consts->tln_mrna_y;
    array[9] = consts->ribo_on_mrna_x;
    array[10] = consts->ribo_off_mrna_x;
    array[11] = consts->tln_mrna_x;
    array[12] = consts->gene_x_kd_y;
    array[13] = consts->gene_x_km_rnap;
    array[14] = consts->gene_x_txn;
    array[15] = consts->gene_y_kd_y;
    array[16] = consts->gene_y_kd_x;
    array[17] = consts->gene_y_basal_km_rnap;
    array[18] = consts->gene_y_basal_txn;
    array[19] = consts->gene_y_active_km_rnap;
    array[20] = consts->gene_y_active_txn;

}

typedef struct OptimOpts{
    
    PatternSearchOpts patopts;
    double upper_bound;
    double lower_bound;
    double wavelength;
    double max_expr;
    double min_expr;
    double burn_in;
    double tend;
    double tstep;
    int report_iters;
    Species* init_values;
    Constants* min;
    Constants* max;

} OptimOpts;

OptimOpts defaultOptimOpts( void ) {
    OptimOpts opts = {
        DEFAULT_PATTERN_SEARCH_OPTS,
        -1,
        -1,
        300,
        6000,
        1000,
        150,
        DEFAULT_TEND,
        DEFAULT_TSTEP,
        DEFAULT_REPORT_ITERS,
        &initValues,
        NULL,
        NULL
    };
    opts.patopts.report = stderr;
    return opts;
};

inline double target( double time, double burn_in, double min_expr, double amp, double wavelength ){
    
    //start time after burn in ends
    return amp/2 - amp*cos( 2*M_PI*( time - burn_in ) / wavelength ) / 2 + min_expr;

}

inline double evaluate( Simulation* sim, OptimOpts* opts ){

    double burn_in = opts->burn_in;
    double min_expr = opts->min_expr;
    double amp = opts->max_expr - min_expr;
    double wavelength = opts->wavelength;

    size_t n = sim->n;
    Species* vals = sim->vals;
    double* time = sim->time;

    double err = 0;
    size_t timepoints = 0;
    size_t i;
    for( i = 0 ; i < n ; ++i ){
        if( time[i] < burn_in ){ //DON'T PENALIZE BURN IN PERIOD
            continue;
        }
        double expect = target( time[i], burn_in, min_expr, amp, wavelength );
        double actual = vals[i].y;
        err += fabs( expect - actual );
        ++timepoints;
    }

    return err / timepoints;

}

double objective( const double* params, size_t n, void* optimopts ){

    OptimOpts* opts = ( OptimOpts* ) optimopts;

    Constants consts;
    arrayToConstants( params, &consts );
    Simulation sim = simulate( &consts, opts->init_values, opts->tend, opts->tstep, NULL, opts->report_iters );
    double err = evaluate( &sim, opts );
    free( sim.vals );
    free( sim.time );

    return err;

}

void optimize( Constants* consts, Species* init_vals, OptimOpts* opts ){
    
    double v0[21];
    constantsToArray( consts, v0 );
    double result[21];
    
    double ub[21];
    if( opts->max != NULL ){
        
        constantsToArray( opts->max, ub );
        opts->patopts.upper_bound = ub;

    }else if( opts->upper_bound >= 0 ){

        double mult = opts->upper_bound;
        size_t i;
        for( i = 0 ; i < 21 ; ++i ){
            ub[i] = mult*v0[i];
        }
        opts->patopts.upper_bound = ub;

    }

    double lb[21];
    if( opts->min != NULL ){

        constantsToArray( opts->min, lb );
        opts->patopts.lower_bound = lb;

    }else if( opts->lower_bound >= 0 ){

        double mult = opts->lower_bound;
        size_t i;
        for( i = 0 ; i < 21 ; ++i ){
            lb[i] = mult * v0[i];
        }
        opts->patopts.lower_bound = lb;

    }

    double best = patternSearch( &objective, v0, result, 21, opts, &(opts->patopts) );
    
    arrayToConstants( result, consts );

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

typedef enum Mode { SIMULATE, OPTIMIZE } Mode;

int assignConstants( Constants* consts, const char* tag, double value ){
    //massive if else for checking the tag
    if( strcmp( tag, "deg_x" ) == 0 ){
        consts->deg_x = value;
    }else if( strcmp( tag, "deg_y" ) == 0 ){
        consts->deg_y = value;
    }else if( strcmp( tag, "coop_y" ) == 0 ){
        consts->coop_y = value;
    }else if( strcmp( tag, "coop_x" ) == 0 ){
        consts->coop_x = value;
    }else if( strcmp( tag, "deg_mrna_y" ) == 0 ){
        consts->deg_mrna_y = value;
    }else if( strcmp( tag, "deg_mrna_x" ) == 0 ){
        consts->deg_mrna_x = value;
    }else if( strcmp( tag, "ribo_on_mrna_y" ) == 0 ){
        consts->ribo_on_mrna_y = value;
    }else if( strcmp( tag, "ribo_off_mrna_y" ) == 0 ){
        consts->ribo_off_mrna_y = value;
    }else if( strcmp( tag, "tln_mrna_y" ) == 0 ){
        consts->tln_mrna_y = value;
    }else if( strcmp( tag, "ribo_on_mrna_x" ) == 0 ){
        consts->ribo_on_mrna_x = value;
    }else if( strcmp( tag, "ribo_off_mrna_x" ) == 0 ){
        consts->ribo_off_mrna_x = value;
    }else if( strcmp( tag, "tln_mrna_x" ) == 0 ){
        consts->tln_mrna_x = value;
    }else if( strcmp( tag, "gene_x_kd_y" ) == 0 ){
        consts->gene_x_kd_y = value;
    }else if( strcmp( tag, "gene_x_km_rnap" ) == 0 ){
        consts->gene_x_km_rnap = value;
    }else if( strcmp( tag, "gene_x_txn" ) == 0 ){
        consts->gene_x_txn = value;
    }else if( strcmp( tag, "gene_y_kd_y" ) == 0 ){
        consts->gene_y_kd_y = value;
    }else if( strcmp( tag, "gene_y_kd_x" ) == 0 ){
        consts->gene_y_kd_x = value;
    }else if( strcmp( tag, "gene_y_basal_km_rnap" ) == 0 ){
        consts->gene_y_basal_km_rnap = value;
    }else if( strcmp( tag, "gene_y_basal_txn" ) == 0 ){
        consts->gene_y_basal_txn = value;
    }else if( strcmp( tag, "gene_y_active_km_rnap" ) == 0 ){
        consts->gene_y_active_km_rnap = value;
    }else if( strcmp( tag, "gene_y_active_txn" ) == 0 ){
        consts->gene_y_active_txn = value;
    }else{
        //unrecognized tag return failure
        fprintf( stderr, "Unrecognized constants tag: %s\n", tag );
        return 0;
    }
    return 1;
}

int parseConstants( FILE* f, Constants* consts ){
    
    int n;
    char tag[20];
    double val;
    while( ( n = fscanf( f, " %19s = %lf\n", tag, &val ) ) == 2 ){
        if( !assignConstants( consts, tag, val ) ){
            return 0;
        }
    }
    return 1;

}

void printConstants( FILE* f, Constants* consts ){

    fprintf( f, "%s = %lf\n", "deg_y", consts->deg_y);
    fprintf( f, "%s = %lf\n", "deg_x", consts->deg_x);
    fprintf( f, "%s = %lf\n", "coop_y", consts->coop_y);
    fprintf( f, "%s = %lf\n", "coop_x", consts->coop_x);

    fprintf( f, "%s = %lf\n", "deg_mrna_y", consts->deg_mrna_y);
    fprintf( f, "%s = %lf\n", "deg_mrna_x", consts->deg_mrna_x);

    fprintf( f, "%s = %lf\n", "ribo_on_mrna_y", consts->ribo_on_mrna_y);
    fprintf( f, "%s = %lf\n", "ribo_off_mrna_y", consts->ribo_off_mrna_y);
    fprintf( f, "%s = %lf\n", "tln_mrna_y", consts->tln_mrna_y);

    fprintf( f, "%s = %lf\n", "ribo_on_mrna_x", consts->ribo_on_mrna_x);
    fprintf( f, "%s = %lf\n", "ribo_off_mrna_x", consts->ribo_off_mrna_x);
    fprintf( f, "%s = %lf\n", "tln_mrna_x", consts->tln_mrna_x);

    fprintf( f, "%s = %lf\n", "gene_x_kd_y", consts->gene_x_kd_y);
    fprintf( f, "%s = %lf\n", "gene_x_km_rnap", consts->gene_x_km_rnap);
    fprintf( f, "%s = %lf\n", "gene_x_txn", consts->gene_x_txn);

    fprintf( f, "%s = %lf\n", "gene_y_kd_y", consts->gene_y_kd_y);
    fprintf( f, "%s = %lf\n", "gene_y_kd_x", consts->gene_y_kd_x);
    fprintf( f, "%s = %lf\n", "gene_y_basal_km_rnap", consts->gene_y_basal_km_rnap);
    fprintf( f, "%s = %lf\n", "gene_y_basal_txn", consts->gene_y_basal_txn);
    fprintf( f, "%s = %lf\n", "gene_y_active_km_rnap", consts->gene_y_active_km_rnap);
    fprintf( f, "%s = %lf\n", "gene_y_active_txn", consts->gene_y_active_txn);

}

int parseSpeciesVals( FILE* f, Species* spec ){

    int n;
    char tag[20];
    double val;
    while( ( n = fscanf( f, " %19s = %lf\n", tag, &val ) ) == 2 ){
        if( strcmp( tag, "rnap" ) == 0 ){
            spec->rnap = val;
        }
        else if( strcmp( tag, "ribo" ) == 0 ){
            spec->ribo = val;
        }
        else if( strcmp( tag, "y" ) == 0 ){
            spec->y = val;
        }
        else if( strcmp( tag, "x" ) == 0 ){
            spec->x = val;
        }
        else if( strcmp( tag, "mrna_y" ) == 0 ){
            spec->mrna_y = val;
        }
        else if( strcmp( tag, "mrna_x" ) == 0 ){
            spec->mrna_x = val;
        }
        else if( strcmp( tag, "gene_y" ) == 0 ){
            spec->gene_y = val;
        }
        else if( strcmp( tag, "gene_x" ) == 0 ){
            spec->gene_x = val;
        }else{
            fprintf( stderr, "Unrecognized initial values tag: %s\n", tag );
            return 0;
        }
    }

    //success
    return 1;

}

int parseOptimOpts( FILE* f, OptimOpts* opts ){

    int n;
    char tag[20];
    double val;
    while( (n = fscanf( f, " %19s = %lf\n", tag, &val ) ) == 2 ){
        //switch on tag to fill values
        if( strcmp( tag, "upper_bound" ) == 0 ){
            opts->upper_bound = val;
        }
        else if( strcmp( tag, "lower_bound" ) == 0 ){
            opts->lower_bound = val;
        }
        else if( strcmp( tag, "mesh_tol" ) == 0 ){
            opts->patopts.mesh_tol = val;
        }
        else if( strcmp( tag, "fnc_tol" ) == 0 ){
            opts->patopts.fnc_tol = val;
        }
        else if( strcmp( tag, "mesh_scale" ) == 0 ){
            opts->patopts.mesh_scale = val;
        }
        else if( strcmp( tag, "mesh_size_init" ) == 0 ){
            opts->patopts.mesh_size_init = val;
        }
        else if( strcmp( tag, "max_iters" ) == 0 ){
            opts->patopts.max_iters = val;
        }
        else if( strcmp( tag, "wavelength" ) == 0 ){
            opts->wavelength = val;
        }
        else if( strcmp( tag, "max_expr" ) == 0 ){
            opts->max_expr = val;
        }
        else if( strcmp( tag, "min_expr" ) == 0 ){
            opts->min_expr = val;
        }
        else if( strcmp( tag, "burn_in" ) == 0 ){
            opts->burn_in = val;
        }else{
            fprintf( stderr, "Unrecognized tag: %s\n", tag );
            return 0;
        }
    }

    return 1;

}

void usage(){
    
    fprintf( stderr, "Usage: simulate [-h] [-p FILE] [-i FILE] [-o FILE] [-t FLOAT] [-s FLOAT] [-r INT] [-m STRING] [--optimopts FILE] [--min FILE] [--max FILE]\n" );
    fprintf( stderr, "\n-p\t--parameters\tFile containing model parameters\n" );
    fprintf( stderr, "-i\t--initialValues\tFile containing initial values of model components\n" );
    fprintf( stderr, "-t\t--timespan\tTimespan to simulate, simulation runs from t=0 to t=timespan, default=1000\n" );
    fprintf( stderr, "-s\t--timestep\tStep size of the simulation, time t_i = t_{i-1} - timestep, default=0.01\n");
    fprintf( stderr, "-r\t--report\tComponent values are output every report iterations of the simulation, default=100\n" );
    fprintf( stderr, "-m\t--mode\tMode to run, options are `simulate' and `optimize', default=`simulate'\n" );
    fprintf( stderr, "\t--optimopts\tFile specifying optimization options\n" );
    fprintf( stderr, "\t--min\tFile specifying minimum value of constants over which to optimize\n" );
    fprintf( stderr, "\t--max\tFile specifying maximum value of constants over which to optimize\n" );
    fprintf( stderr, "-o\t--output\tFile to output the simulation results to, default=stdout\n" );
    fprintf( stderr, "-h\t--help\tDisplays this help message\n" );

}


int main( int argc, char* argv[] ){

    FILE* out = stdout;
    char* param_file = NULL;
    char* initvals_file = NULL;
    char* out_file = NULL;
    char* optimopts = NULL;
    char* min_file = NULL;
    char* max_file = NULL;
    double tend = DEFAULT_TEND;
    double tstep = DEFAULT_TSTEP;
    int report_iters = DEFAULT_REPORT_ITERS;
    Mode mode = SIMULATE;

    static struct option param_opts[] = 
    {
        { "parameters", required_argument, 0, 'p' },
        { "initialValues", required_argument, 0, 'i' },
        { "output", required_argument, 0, 'o' },
        { "timespan", required_argument, 0, 't' },
        { "timestep", required_argument, 0, 's' },
        { "report", required_argument, 0, 'r' },
        { "mode", required_argument, 0, 'm' },
        { "optimopts", required_argument, 0, '1' },
        { "min", required_argument, 0, '2' },
        { "max", required_argument, 0, '3' },
        { "help", no_argument, 0, 'h' },
        { 0, 0, 0, 0 }
    };

    //parse command line arguments
    int opt_idx = 0;
    int c;
    while( ( c = getopt_long( argc, argv, "p:i:o:t:s:r:m:h", 
                    param_opts, &opt_idx ) ) != -1 ){

        switch( c ){
            case 'p':
                param_file = optarg;
                break;
            case 'i':
                initvals_file = optarg;
                break;
            case 'o':
                out_file = optarg; 
                break;
            case 't':
                tend = atof( optarg );
                break;
            case 's':
                tstep = atof( optarg );
                break;
            case 'r':
                report_iters = atoi( optarg );
                break;
            case 'm':
                if( strcmp( optarg, "simulate" ) == 0 ){
                    mode = SIMULATE;
                }else if( strcmp( optarg, "optimize" ) == 0 ){
                    mode = OPTIMIZE;
                }else{
                    fprintf( stderr, "Unrecognized mode: %s\n", optarg );
                    fprintf( stderr,
                        "Valid modes are `%s', `%s'\n",
                        "simulate",
                        "optimize" );
                    return 1;
                }
                break;
            case 'h':
                usage();
                return 0;
            case '1':
                optimopts = optarg;
                break;
            case '2':
                min_file = optarg;
                break;
            case '3':
                max_file = optarg;
                break;
        }

    }

    //validate command line arguments

    size_t num;
    Constants params = initParams;
    Species init_vals = initValues;
    OptimOpts optimparams = defaultOptimOpts( );
    optimparams.tend = tend;
    optimparams.tstep = tstep;
    optimparams.report_iters = report_iters;
    optimparams.init_values = &init_vals;
    Constants min;
    Constants max;

    if( param_file != NULL ){
        FILE* f = fopen( param_file, "r" );
        if( f == NULL ){
            fprintf( stderr, "Unable to open parameters file: %s\n", param_file );
            return 1;
        }
        if( !parseConstants( f, &params ) ){
            //error occurred parsing the constants file
            fprintf( stderr, "Error parsing parameters file: %s\n", param_file );
            return 1;
        }
        fclose( f );
    }

    if( initvals_file != NULL ){
        FILE* f = fopen( initvals_file, "r" );
        if( f == NULL ){
            fprintf( stderr, "Unable to open initial values file: %s\n", initvals_file );
            return 1;
        }
        if( !parseSpeciesVals( f, &init_vals ) ){
            return 1;
        }
        fclose( f );
    }


    if( out_file != NULL ){
        out = fopen( out_file, "w" );
        if( out == NULL ){
            fprintf( stderr, "Unable to open output file: %s\n", optarg );
            return 1;
        }
    }

    if( tend < 0 ){
        fprintf( stderr, "Timespan must be >= 0 but was %lf\n", tend );
        return 1;
    }

    if( tstep >= tend ){
        fprintf( stderr, "Timestep must be < timespan but was %lf\n", tstep );
        return 1;
    }

    if( optimopts != NULL ){
        FILE* f = fopen( optimopts, "r" );
        if( f == NULL ){
            fprintf( stderr, "Unable to open optimization options file: %s\n", optimopts );
            return 1;
        }
        if( !parseOptimOpts( f, &optimparams ) ){
            return 1;
        }
        fclose( f );
    }

    if( min_file != NULL ){
        FILE* f = fopen( min_file, "r" );
        if( f == NULL ){
            fprintf( stderr, "Unable to open minimum parameters file: %s\n", min_file );
            return 1;
        }
        if( !parseConstants( f, &min ) ){
            return 1;
        }
        optimparams.min = &min;
        fclose( f );
    }

    if( max_file != NULL ){
        FILE* f = fopen( max_file, "r" );
        if( f == NULL ){
            fprintf( stderr, "Unable to open maximum parameters file: %s\n", max_file );
            return 1;
        }
        if( !parseConstants( f, &max ) ){
            return 1;
        }
        optimparams.max = &max;
        fclose( f );
    }
    
   
    fprintf( stderr, "Running simulation with TEND = %lf, TSTEP = %lf, REPORT_ITERS = %d\n", tend, tstep, report_iters
           );

    switch( mode ){
        case SIMULATE:
            simulate( &params, &init_vals, tend, tstep, out, report_iters );
            break;
        case OPTIMIZE:
            optimize( &params, &init_vals, &optimparams );
            Simulation sim = simulate( &params, &init_vals, tend, tstep, NULL, report_iters );
            printConstants( stderr, &params );
            fprintf( out, "Time X Y mRNA_X mRNA_Y Objective\n" );
            size_t i;
            double burn_in = optimparams.burn_in;
            double min_expr = optimparams.min_expr;
            double amp = optimparams.max_expr - min_expr;
            double wavelength = optimparams.wavelength;
            for( i = 0 ; i < sim.n ; ++i ){
                Species cur = sim.vals[i];
                double time = sim.time[i];
                double obj = target( time, burn_in, min_expr, amp, wavelength );
                fprintf( out, "%lf %lf %lf %lf %lf %lf\n", time, cur.x, cur.y, cur.mrna_x, cur.mrna_y, obj );
            }
            free( sim.vals );
            free( sim.time );
            break;
    }

    return 0;


}

