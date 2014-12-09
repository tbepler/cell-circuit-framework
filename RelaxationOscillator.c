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

    //ensure values never less than zero
    if( next->y < 0 ){
        next->y = 0;
    }
    if( next->x < 0 ){
        next->x = 0;
    }

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
    //ensure values never less than zero
    if( next->mrna_x < 0 ){
        next->mrna_x = 0;
    }
    if( next->mrna_y < 0 ){
        next->mrna_y = 0;
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
    size_t n = tend / tstep / reportIters + 1;
    sim.vals = malloc( n * sizeof( Species ) );
    sim.time = malloc( n * sizeof( double ) );
    sim.n = n;
    
    memcpy( sim.vals, &cur, sizeof( *sim.vals ) );
    sim.time[0] = time;
    if( out != NULL ){
        report( out, &cur, time );
    }
    
    size_t idx = 1;

    while( time < tend ){
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

#define WL 300
#define AMP 8000

inline double target( double time ){
    
    return AMP/2 + AMP*cos( 2*M_PI*time / WL ) / 2;

}

inline double evaluate( Simulation* sim ){

    size_t n = sim->n;
    Species* vals = sim->vals;
    double* time = sim->time;
    double err = 0;
    size_t i;
    for( i = 0 ; i < n ; ++i ){
        double expect = target( time[i] );
        double actual = vals[i].y;
        err += fabs( expect - actual );
    }

    return err / n;

}

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

double objective( const double* params, size_t n ){

    Constants consts;
    arrayToConstants( params, &consts );
    Simulation sim = simulate( &consts, &initValues, 1000, 0.1, NULL, 10 );
    double err = evaluate( &sim );
    free( sim.vals );
    free( sim.time );

    return err;

}

void optimize( Constants* consts ){
    
    double v0[21];
    constantsToArray( consts, v0 );
    double result[21];
    double ub[21];
    size_t i;
    for( i = 0 ; i < 21 ; ++i ){
        ub[i] = 10*v0[i];
    }
    double lb[21] = {0};
    double best = patternSearch( &objective, v0, result, 21, ub, lb, stderr );
    
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
        return 0;
    }
    return 1;
}

char* nextLine( FILE* file ){

    char* line;
    while( ( line = readLine( file ) ) != NULL &&
            ( line[0] == '\0' || line[0] == '#' ) ){
        // # marks comment line and ignore empty lines
        free( line ); //need to free string malloc'd by readLine
    }
    return line; 

}

int parseConstants( FILE* file, Constants* consts, char** tags, size_t n_tags ){

    char* line = nextLine( file );
    if( line == NULL ){
        //EOF reached but no entries
        return 0;
    }

    //convert line into vector of doubles
    size_t n;
    double* vec = parseDoubles( line, DELIMS, &n );
    //free the line string since no longer needed
    free( line );
    if( n < n_tags ){
        //error, not enough elements on this line to specify all tags
        free( vec );
        return 0;
    }

    size_t i;
    for( i = 0 ; i < n ; ++i ){
        //assign fields of consts
        if( !assignConstants( consts, tags[i], vec[i] ) ){
            //error occurred on assignment
            free( vec );
            return 0;
        }
    }
    
    //all assignments completed
    free( vec );
    return 1;

}

char** parseHeader( FILE* file, size_t* p_n ){

    //read first nonempty non # commented line
    char* line = nextLine( file );
    if( line == NULL ){
        //EOF reached but no entries
        return NULL;
    } 
    
    char** header = split( line, DELIMS, p_n );
    free( line );
    return header;

}

void printStrArray( FILE* out, char** strs, size_t p_n ){
    
    size_t i;
    for( i = 0 ; i < p_n - 1 ; ++i ){
        fprintf( out, "%s ", strs[i] );
    }
    fprintf( out, "%s\n", strs[i] );

}

void printConstants( FILE* out, const Constants* consts ){
    
    fprintf( out, "%lf ", consts->deg_y );
    fprintf( out, "%lf ", consts->deg_x );
    fprintf( out, "%lf ", consts->coop_y );
    fprintf( out, "%lf ", consts->coop_x );

    fprintf( out, "%lf ", consts->deg_mrna_y );
    fprintf( out, "%lf ", consts->deg_mrna_x );

    fprintf( out, "%lf ", consts->ribo_on_mrna_y );
    fprintf( out, "%lf ", consts->ribo_off_mrna_y );
    fprintf( out, "%lf ", consts->tln_mrna_y );

    fprintf( out, "%lf ", consts->ribo_on_mrna_x );
    fprintf( out, "%lf ", consts->ribo_off_mrna_x );
    fprintf( out, "%lf ", consts->tln_mrna_x );

    fprintf( out, "%lf ", consts->gene_x_kd_y );
    fprintf( out, "%lf ", consts->gene_x_km_rnap );
    fprintf( out, "%lf ", consts->gene_x_txn );

    fprintf( out, "%lf ", consts->gene_y_kd_y );
    fprintf( out, "%lf ", consts->gene_y_kd_x );
    fprintf( out, "%lf ", consts->gene_y_basal_km_rnap );
    fprintf( out, "%lf ", consts->gene_y_basal_txn );
    fprintf( out, "%lf ", consts->gene_y_active_km_rnap );
    fprintf( out, "%lf\n", consts->gene_y_active_txn );

}

void usage(){
    
    fprintf( stderr, "Usage: simulate [-h] [-p FILE] [-i FILE] [-o FILE] [-t FLOAT] [-s FLOAT] [-r INT] [-m STRING]\n" );
    fprintf( stderr, "\n-p\t--parameters\tFile containing model parameters\n" );
    fprintf( stderr, "-i\t--initialValues\tFile containing initial values of model components\n" );
    fprintf( stderr, "-t\t--timespan\tTimespan to simulate, simulation runs from t=0 to t=timespan, default=1000\n" );
    fprintf( stderr, "-s\t--timestep\tStep size of the simulation, time t_i = t_{i-1} - timestep, default=0.01\n");
    fprintf( stderr, "-r\t--report\tComponent values are output every report iterations of the simulation, default=100\n" );
    fprintf( stderr, "-m\t--mode\tMode to run, options are `simulate' and `optimize', default=`simulate'\n" );
    fprintf( stderr, "-o\t--output\tFile to output the simulation results to, default=stdout\n" );
    fprintf( stderr, "-h\t--help\tDisplays this help message\n" );

}


int main( int argc, char* argv[] ){

    FILE* out = stdout;
    char* param_file = NULL;
    char* initvals_file = NULL;
    char* out_file = NULL;
    double tend = 1000;
    double tstep = 0.01;
    int report_iters = 100;
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
        }

    }

    //validate command line arguments

    size_t num;
    char** header;
    Constants params = initParams;
    Species init_vals = initValues;

    if( param_file != NULL ){
        FILE* f = fopen( param_file, "r" );
        if( f == NULL ){
            fprintf( stderr, "Unable to open parameters file: %s\n", param_file );
            return 1;
        }
        header = parseHeader( f, &num );
        if( header == NULL ){
            //error occurred
            fprintf( stderr, "Error parsing parameters file: %s\n", param_file );
            return 1;
        }
        if( !parseConstants( f, &params, header, num ) ){
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
        //TODO process file
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
    
   
    fprintf( stderr, "Running simulation with TEND = %lf, TSTEP = %lf, REPORT_ITERS = %d\n", tend, tstep, report_iters
           );

    switch( mode ){
        case SIMULATE:
            simulate( &params, &init_vals, tend, tstep, out, report_iters );
            break;
        case OPTIMIZE:
            optimize( &params );
            simulate( &params, &init_vals, tend, tstep, out, report_iters );
            break;
    }

    return 0;


}

