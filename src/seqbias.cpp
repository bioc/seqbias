
#include "sequencing_bias.hpp"
#include <algorithm>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <cstdio>

extern "C" {

#include <htslib/faidx.h>
#include <samtools-1.7-compat.h>


#ifdef ERROR
#undef ERROR
#endif


// more R namespace pollution
#if defined(nrows)
#undef nrows
#endif

#if defined(ncols)
#undef ncols
#endif


typedef struct {
    samfile_t*   f;
    bam_index_t* idx;
} indexed_bam_f;


/* several functions need to covert arguments that give genomic coordinates */
void coerce_genomic_coords( SEXP seqname,
                            SEXP start,
                            SEXP end,
                            SEXP strand,
                            const char** c_seqname,
                            pos_t*       c_start,
                            pos_t*       c_end,
                            strand_t*    c_strand )
{
    if( !IS_CHARACTER(seqname) || LENGTH(seqname) != 1 ) error( "seqname must be character(1)" );
    *c_seqname = translateChar( STRING_ELT( seqname, 0 ) );

    if( !IS_INTEGER(start) ) error( "start is non-integer" );
    *c_start = asInteger(start);
    if( *c_start < 0 ) error( "start must be positive" );

    if( !IS_INTEGER(end) ) error( "end is non-integer" );
    *c_end = asInteger(end);
    if( *c_end < 0 ) error( "end must be positive" );

    if( *c_end < *c_start ) error( "end must be greater or equal to start" );

    if( !IS_CHARACTER(strand) || LENGTH(strand) != 1 ) error( "strand must be character(1)" );
    const char* c_strand_str = translateChar( STRING_ELT( strand, 0 ) );

    if( strlen(c_strand_str) != 1 ) error( "strand should be be one character" );

    if( c_strand_str[0] == '+' )      *c_strand = strand_pos;
    else if( c_strand_str[0] == '-' ) *c_strand = strand_neg;
    else                              *c_strand = strand_na;
}


static void free_seqbias( SEXP R_seqbias )
{
    delete (sequencing_bias*)EXTPTR_PTR(R_seqbias);
}

SEXP seqbias_load( SEXP ref_fn, SEXP model_fn )
{
    if( !IS_CHARACTER(ref_fn) || LENGTH(ref_fn) != 1 ) error( "'ref_fn' must be character(1)" );
    if( !IS_CHARACTER(model_fn) || LENGTH(model_fn) != 1 ) error( "'model_fn' must be character(1)" );

    const char* c_ref_fn   = translateChar( STRING_ELT( ref_fn, 0 ) );
    const char* c_model_fn = translateChar( STRING_ELT( model_fn, 0 ) );

    sequencing_bias* c_seqbias = new sequencing_bias( c_ref_fn, c_model_fn );

    SEXP seqbias = R_MakeExternalPtr( (void*)c_seqbias, R_NilValue, R_NilValue );
    R_RegisterCFinalizer( seqbias, free_seqbias );

    return seqbias;
}


SEXP seqbias_save( SEXP seqbias, SEXP fn )
{
    sequencing_bias* c_seqbias;
    if( TYPEOF(seqbias) != EXTPTRSXP || !(c_seqbias = (sequencing_bias*)EXTPTR_PTR(seqbias)) ) {
        error( "first argument is not a proper seqbias class." );
    }

    if( !IS_CHARACTER(fn) || LENGTH(fn) != 1 ) error( "'fn' must be character(1)" );
    const char* c_fn = translateChar( STRING_ELT( fn, 0 ) );

    c_seqbias->save_to_file( c_fn );

    return R_NilValue;
}


SEXP seqbias_fit( SEXP ref_fn,
                  SEXP reads_fn,
                  SEXP n,
                  SEXP L, SEXP R )
{
    if( !IS_CHARACTER(ref_fn)   || LENGTH(ref_fn) != 1 )   error( "'ref_fn' must be character(1)" );
    if( !IS_CHARACTER(reads_fn) || LENGTH(reads_fn) != 1 ) error( "'reads_fn' must be character(1)" );
    if( !IS_INTEGER(n) ) error( "'n' is non-integer" );
    if( !IS_INTEGER(L) ) error( "'L' is non-integer" );
    if( !IS_INTEGER(R) ) error( "'R' is non-integer" );

    const char* c_ref_fn   = translateChar( STRING_ELT( ref_fn, 0 ) );
    const char* c_reads_fn = translateChar( STRING_ELT( reads_fn, 0 ) );

    int c_n = asInteger( n );
    int c_L = asInteger( L );
    int c_R = asInteger( R );

    if( c_n < 1 ) error( "'n' must be positive" );
    if( c_L < 0 ) error( "'L' must be non-negative" );
    if( c_R < 0 ) error( "'R' must be non-negative" );

    sequencing_bias* seqbias = new sequencing_bias( c_ref_fn, c_reads_fn, c_n, c_L, c_R );

    SEXP R_seqbias = R_MakeExternalPtr( (void*)seqbias, R_NilValue, R_NilValue );
    R_RegisterCFinalizer( R_seqbias, free_seqbias );

    return R_seqbias;
}


SEXP seqbias_predict( SEXP seqbias,
                      SEXP seqname,
                      SEXP start,
                      SEXP end,
                      SEXP strand )
{
    sequencing_bias* c_seqbias;
    if( TYPEOF(seqbias) != EXTPTRSXP || !(c_seqbias = (sequencing_bias*)EXTPTR_PTR(seqbias)) ) {
        error( "first argument is not a proper seqbias class." );
    }


    const char* c_seqname;
    pos_t c_start, c_end;
    strand_t c_strand;


    coerce_genomic_coords( seqname, start, end, strand,
                           &c_seqname, &c_start, &c_end, &c_strand );

    /* we expect 1-based coordinates, but work with 0-based */
    c_start -= 1;
    c_end   -= 1;
                           
    if( c_strand != 0 && c_strand != 1 ) {
        warning( "strand should be '+' or '-'" );
        return R_NilValue;
    }

    size_t n = c_end - c_start + 1;
    double* c_v = c_seqbias->get_bias( c_seqname, c_start, c_end, c_strand );
    SEXP v;
    PROTECT( v = allocVector( REALSXP, c_end - c_start + 1 ) );

    size_t i;
    for( i = 0; i < n; i++ ) REAL(v)[i] = c_v[i];

    delete[] c_v;

    UNPROTECT(1);
    return v;
}



void seqbias_close_bam( SEXP bam_ptr )
{
    if( TYPEOF(bam_ptr) != EXTPTRSXP ) error( "argument is not a indexed bam pointer" );
    indexed_bam_f* c_bam_ptr = (indexed_bam_f*)EXTPTR_PTR( bam_ptr );

    bam_index_destroy( c_bam_ptr->idx );
    samclose( c_bam_ptr->f );

    free( c_bam_ptr );
    c_bam_ptr = NULL;
}


SEXP seqbias_open_bam( SEXP reads_fn )
{
    if( !IS_CHARACTER(reads_fn) || LENGTH(reads_fn) != 1 ) {
        error( "'reads_fn' must be character(1)" );
    }
    const char* c_reads_fn = translateChar( STRING_ELT( reads_fn, 0 ) );

    samfile_t* f;
    bam_index_t* idx;

    f = samopen( c_reads_fn, "rb", NULL ); 
    if( f == NULL ) error( "can't open BAM file" );

    idx = bam_index_load( c_reads_fn );
    if( idx == NULL ) error( "can't open BAM index file" );

    indexed_bam_f* ib = (indexed_bam_f*)malloc( sizeof(indexed_bam_f) );
    ib->f   = f;
    ib->idx = idx;

    SEXP ib_ptr = R_MakeExternalPtr( (void*)ib, R_NilValue, R_NilValue );
    R_RegisterCFinalizer( ib_ptr, seqbias_close_bam );

    return ib_ptr;
}


SEXP seqbias_count_reads( SEXP seqbias,
                          SEXP bam_ptr,
                          SEXP seqname,
                          SEXP start,
                          SEXP end,
                          SEXP strand,
                          SEXP sum_counts )
{
    if( TYPEOF(bam_ptr) != EXTPTRSXP ) error( "argument is not a indexed bam pointer" );
    indexed_bam_f* c_bam_ptr = (indexed_bam_f*)EXTPTR_PTR( bam_ptr );

    const char* c_seqname;
    pos_t c_start, c_end;
    strand_t c_strand;

    coerce_genomic_coords( seqname, start, end, strand,
                           &c_seqname, &c_start, &c_end, &c_strand );


    /* we expect 1-based coordinates, but use 0-based */
    c_start -= 1;
    c_end   -= 1;
                           

    double* bs[2] = {NULL, NULL};

    if( !isNull(seqbias) ) {

        SEXP idx;
        PROTECT( idx = allocVector( STRSXP, 1 ) );
        SET_STRING_ELT(idx, 0, mkChar("ptr"));

        sequencing_bias* c_seqbias = NULL;
        if( TYPEOF(GET_SLOT(seqbias, idx)) != EXTPTRSXP ||
            !(c_seqbias = (sequencing_bias*)EXTPTR_PTR(GET_SLOT(seqbias, idx))) ) {
            error( "first argument is not a proper seqbias class." );
        }
        
        if (c_strand == strand_na || c_strand == strand_pos ) {
            bs[0] = c_seqbias->get_bias(c_seqname, c_start, c_end, strand_pos);
        }

        if (c_strand == strand_na || c_strand == strand_neg ) {
            bs[1] = c_seqbias->get_bias(c_seqname, c_start, c_end, strand_neg);
            std::reverse(bs[1], bs[1] + (c_end - c_start + 1));
        }

        UNPROTECT(1);
    }

    bool c_sum_counts = asLogical(sum_counts) == TRUE;

    /* init vector */
    SEXP v;
    if (c_sum_counts) {
        PROTECT( v = allocVector( REALSXP, 1 ) );
        REAL(v)[0] = 0;
    }
    else {
        PROTECT( v = allocVector( REALSXP, c_end - c_start + 1 ) );
        pos_t i;
        for( i = 0; i < c_end - c_start + 1; i++ ) REAL(v)[i] = 0;
    }


    const size_t region_len = 1024;
    char* region = new char[region_len];

    int bam_ref_id, bam_start, bam_end, err;
    err = snprintf(region, region_len, "%s:%ld-%ld", c_seqname, c_start, c_end );
    err = bam_parse_region( c_bam_ptr->f->header, region,
                            &bam_ref_id, &bam_start, &bam_end );

    delete [] region;

    /* if the region is not present in the bam file index, just return 0's */
    if( err != 0 || bam_ref_id < 0 ) {
        UNPROTECT(1);
        return v;
    }


    bam_iter_t it = bam_iter_query( c_bam_ptr->idx, bam_ref_id, bam_start, bam_end );
    bam1_t* b = bam_init1();
    pos_t x;
    strand_t s;


    while( bam_iter_read( c_bam_ptr->f->x.bam, it, b ) > 0 ) {
        s = (strand_t) bam1_strand(b);
        if( c_strand != strand_na && s !=  c_strand ) continue;

        x = bam1_strand(b) == 1 ? bam_calend( &b->core,  bam1_cigar(b) ) - 1 : b->core.pos;
        if (x < c_start || x > c_end) continue;

        if (c_sum_counts) {
            REAL(v)[0] += bs[s] == NULL ? 1.0 : (1.0 / bs[s][x - c_start]);
        }
        else {
            REAL(v)[x - c_start] += bs[s] == NULL ? 1.0 : (1.0 / bs[s][x - c_start]);
        }
    }

    if( c_strand == strand_neg && !c_sum_counts ) {
        std::reverse(REAL(v), REAL(v) + (c_end - c_start + 1));
    }

    bam_iter_destroy( it );
    bam_destroy1(b);

    delete [] bs[0];
    delete [] bs[1];

    UNPROTECT(1);
    return v;
}




void dealloc_kmer_matrix( SEXP M )
{
    if( TYPEOF(M) != EXTPTRSXP ) error( "argument is not a kmer_matrix pointer" );
    kmer_matrix* c_M = (kmer_matrix*)EXTPTR_PTR(M);
    delete c_M;
}



SEXP seqbias_alloc_kmer_matrix( SEXP n, SEXP k )
{
    int c_n = 0, c_k = 0;

    if( !IS_INTEGER(n) || (c_n = asInteger(n)) <= 0) {
        error( "'n' must be a positive integer" );
    }

    if( !IS_INTEGER(k) || (c_k = asInteger(k)) <= 0) {
        error( "'k' must be a positive integer" );
    }

    kmer_matrix* c_M = new kmer_matrix( c_n, c_k );
    c_M->set_all( 0.0 );

    SEXP M = R_MakeExternalPtr( (void*)c_M, R_NilValue, R_NilValue );
    R_RegisterCFinalizer( M, dealloc_kmer_matrix );

    return M;

}



SEXP seqbias_tally_kmers( SEXP M, SEXP seq, SEXP count, SEXP offset )
{
    if( TYPEOF(M) != EXTPTRSXP ) error( "M is not a kmer_matrix pointer" );
    kmer_matrix* c_M = (kmer_matrix*)EXTPTR_PTR(M);

    if( !IS_CHARACTER(seq) || LENGTH(seq) != 1 ) error( "seq must be character(1)" );
    const char* c_seq = translateChar( STRING_ELT( seq, 0 ) );

    if( !IS_NUMERIC(count) ) error( "count must be numeric" );

    if( !IS_INTEGER(offset) ) error( "offset must be an integer" );
    pos_t c_offset = asInteger(offset);

    size_t n = strlen(c_seq);
    if( (size_t)LENGTH(count) != n ) error( "sequence length mismatches count length" );

    size_t k = c_M->ksize();

    /*
     * Convert the sequence to an array of kmers.
     */

    size_t i;
    kmer kmer_mask = 0;
    for( i = 0; i < k; i++ ) kmer_mask = (kmer_mask<<2) | 0x3;

    kmer* ks = new kmer[n - (k - 1)];
    memset( ks, 0, (n - (k - 1)) * sizeof(kmer) );
    kmer K = 0;
    for( i = 0; i < n; i++ ) {
        K = ((K << 2) | nuc_to_num(c_seq[i])) & kmer_mask;
        if( i >= k-1 ) ks[i-(k-1)] = K;
    }


    /*
     * Walk through the count array tallying kmers.
     */
    pos_t j;
    for( i = 0; i < n; i++ ) {
        if( (pos_t)i >= c_offset &&
            (pos_t)i - c_offset + c_M->nrows() <= n &&
            REAL(count)[i] > 0.0 )
        {
            for( j = 0; (size_t)j < c_M->nrows(); j++ ) {
                (*c_M)( j, ks[(pos_t)i - c_offset + j] ) += REAL(count)[i];
                //(*c_M)( j, ks[i - c_offset + j] ) += 1;
            }
        }
    }


    delete[] ks;
    return R_NilValue;
}


SEXP seqbias_dataframe_from_kmer_matrix( SEXP M, SEXP offset )
{
    if( TYPEOF(M) != EXTPTRSXP ) error( "M is not a kmer_matrix pointer" );
    kmer_matrix* c_M = (kmer_matrix*)EXTPTR_PTR(M);

    if( !IS_INTEGER(offset) ) error( "offset must be an integer" );
    pos_t c_offset = asInteger(offset);

    /* normalize to get a probability distribution */
    c_M->make_distribution();

    size_t n = c_M->nrows();
    size_t m = c_M->ncols();
    size_t k = c_M->ksize();

    SEXP poss, seqs, freqs;
    PROTECT( poss  = allocVector( REALSXP, n*m ) );
    PROTECT( seqs  = allocVector( STRSXP,  n*m ) );
    PROTECT( freqs = allocVector( REALSXP, n*m ) );

    pos_t i;
    kmer K;

    char* Kstr = new char[k+1];

    for( i = 0; (size_t)i < n; i++ ) {
        for( K = 0; K < m; K++ ) {

            /* set pos */
            REAL(poss)[i*m+K] = i - c_offset;

            /* set seq */
            num_to_nuc(Kstr, K, k);
            SET_STRING_ELT( seqs, i*m+K, mkChar( Kstr ) );

            /* set freq */
            REAL(freqs)[i*m+K] = (*c_M)( i, K );
        }
    }
    delete[] Kstr;

    SEXP result;
    PROTECT( result = allocVector( VECSXP, 3 ) );
    SET_VECTOR_ELT( result, 0, poss );
    SET_VECTOR_ELT( result, 1, seqs );
    SET_VECTOR_ELT( result, 2, freqs );
    UNPROTECT(4);

    return result;
}




void R_init_seqbias( DllInfo* info )
{
    R_CallMethodDef methods[] = {
        { "seqbias_fit",         (DL_FUNC) &seqbias_fit, 5 },
        { "seqbias_predict",     (DL_FUNC) &seqbias_predict, 5 },
        { "seqbias_load",        (DL_FUNC) &seqbias_load, 2 },
        { "seqbias_save",        (DL_FUNC) &seqbias_save, 2 },
        { "seqbias_open_bam",       (DL_FUNC) &seqbias_open_bam, 1 },
        { "seqbias_count_reads",    (DL_FUNC) &seqbias_count_reads, 7 },
        { "seqbias_alloc_kmer_matrix",          (DL_FUNC) &seqbias_alloc_kmer_matrix, 2 },
        { "seqbias_tally_kmers",                (DL_FUNC) &seqbias_tally_kmers, 4 },
        { "sebqias_dataframe_from_kmer_matrix", (DL_FUNC) &seqbias_dataframe_from_kmer_matrix, 2 },
        { NULL, NULL, 0 }
    };

    R_registerRoutines( info, NULL, methods, NULL, NULL );
}


}





