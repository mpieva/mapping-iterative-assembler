#ifndef INCLUDED_MYERS_ALIGN_H
#define INCLUDED_MYERS_ALIGN_H

enum Mode { equal, a_is_prefix, b_is_prefix, either_prefix, overlap, primer_special } ;

#ifdef __cplusplus
extern "C" {
int myers_align( const char *seq_a, const char* seq_b, int maxd, Mode mode = equal,
		char*bt=0, int*xout=0, int*yout=0 ) ;

int myers_alignN( const char *seq_a, const char** seqs_b, int maxd, int* ntargets,
		Mode mode = equal, char*bt=0, int*xout=0, int*yout=0 ) ;
}
#else
int myers_align( const char *seq_a, const char* seq_b, int maxd, enum Mode mode,
		char*bt, int*xout, int*yout ) ;

int myers_alignN( const char *seq_a, const char** seqs_b, int maxd, int* ntargets,
		enum Mode mode, char*bt, int*xout, int*yout ) ;
#endif

inline int char_to_bitmap( char x ) 
{
    switch( x & ~32 )
    {
        case 'A': return 1 ;
        case 'C': return 2 ;
        case 'G': return 4 ;
        case 'T': return 8 ;
        case 'U': return 8 ;

        case 'S': return 6 ;
        case 'W': return 9 ;
        case 'R': return 5 ;
        case 'Y': return 10 ;
        case 'K': return 12 ;
        case 'M': return 3 ;

        case 'B': return 14 ;
        case 'D': return 13 ;
        case 'H': return 11 ;
        case 'V': return 7 ;

        case 'N': return 15 ;
        default: return 0 ;
    }
}

inline char complement_char( char x ) 
{
    switch( x & ~32 )
    {
        case 'A': return 'T' ;
        case 'C': return 'G' ;
        case 'G': return 'C' ;
        case 'T': return 'A' ;
        case 'U': return 'A' ;

        case 'S': return 'S' ;
        case 'W': return 'W' ;
        case 'R': return 'Y' ;
        case 'Y': return 'R' ;
        case 'K': return 'M' ;
        case 'M': return 'K' ;

        case 'B': return 'V' ;
        case 'D': return 'H' ;
        case 'H': return 'D' ;
        case 'V': return 'B' ;

        case 'N': return 'N' ;
        default: return x ;
    }
}
                  
#endif
