#ifndef INCLUDED_MYERS_ALIGN
#define INCLUDED_MYERS_ALIGN

enum myers_align_mode { 
	myers_align_globally,
	myers_align_is_prefix,
	myers_align_has_prefix } ;

//! \brief aligns two sequences in O(nd) time
//! This alignment algorithm following Eugene W. Myers: "An O(ND)
//! Difference Algorithm and Its Variations".
//! Both input sequences are ASCIIZ-encoded with IUPAC ambiguity codes.
//! By definition, if ambiguity codes overlap, that's a match, else a
//! mismatch.  Mismatches and gaps count a unit penalty.  If mode is
//! myers_align_globally, both sequences must align completely.  If mode
//! is myers_align_is_prefix, seq_a must align completely as prefix of
//! seq_b.  If mode is myers_align_has_prefix, seq_b must align
//! completely as prefix of seq_a.  
//!
//! Note that the calculation time is O(nd) where n is the length of the
//! best alignment and d the number of differences in it, memory
//! consumption is O(maxd^2).
//!
//! \param seq_a First input sequence.
//! \param mode How to align (i.e. what gaps to count).
//! \param seq_b Second input sequence.
//! \param maxd Maximum penalty to consider.
//! \param bt_a Space to backtrace seq_a into, must have room for
//!             (strlen(seq_a)+maxd+1) characters.
//! \param bt_b Space to backtrace seq_b into, must have room for
//!             (strlen(seq_b)+maxd+1) characters.
//! \return The actual edit distance or UINT_MAX if the edit distance
//!         would be greater than maxd.
//!
unsigned myers_diff( const char *seq_a, enum myers_align_mode mode, const char* seq_b, int maxd, char *bt_a, char *bt_b ) ;

//! \brief converts an IUPAC ambiguity code to a bitmap
//! Each base is represented by a bit, makes checking for matches
//! easier.
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

inline int compatible( char x, char y ) { return (char_to_bitmap(x) & char_to_bitmap(y)) != 0 ; }

inline int min( int a, int b ) { return a < b ? a : b ; }
inline int max( int a, int b ) { return a < b ? b : a ; }
inline int max3( int a, int b, int c ) { return a < b ? max( b, c ) : max( a, c ) ; }

#endif
