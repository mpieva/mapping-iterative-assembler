#include "myers_align.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

inline int match( char a, char b ) { return (char_to_bitmap(a) & char_to_bitmap(b)) != 0 ; }

// [*blech*, this looks and feels like FORTRAN.]
unsigned myers_diff( const char *seq_a, enum myers_align_mode mode, const char* seq_b, int maxd, char *bt_a, char *bt_b ) 
{
	int len_a = strlen( seq_a ), len_b = strlen( seq_b ) ;
	if( maxd > len_a + len_b ) maxd = len_a + len_b ;

	// in vee[d][k], d runs from 0 to maxd; k runs from -d to +d
	int **vee = calloc( maxd, sizeof(int*) ) ;

	int d, dd, k, x, y, r = UINT_MAX ;
	int *v_d_1 = 0, *v_d = 0 ; 															// "array slice" vee[.][d-1]
	for( d = 0 ; d != maxd ; ++d, v_d_1 = v_d )									// D-paths in order of increasing D
	{
		v_d = d + (vee[d] = malloc( (2 * d + 1) * sizeof( int ) )) ; 		// "array slice" vee[.][d]

		for( k = max(-d,-len_a) ; k <= min(d,len_b) ; ++k ) 					// diagonals
		{
			if( d == 0 )         x = 0 ;
			else if(d==1&&k==0)  x =                       v_d_1[ k ]+1 ;
			else if( k == -d   ) x =                                     v_d_1[ k+1 ] ;
			else if( k ==  d   ) x =       v_d_1[ k-1 ]+1 ;									// argh, need to check for d first, b/c -d+2 could be equal to d
			else if( k == -d+1 ) x = max(                  v_d_1[ k ]+1, v_d_1[ k+1 ] ) ;
			else if( k ==  d-1 ) x = max(  v_d_1[ k-1 ]+1, v_d_1[ k ]+1 ) ;
			else                 x = max3( v_d_1[ k-1 ]+1, v_d_1[ k ]+1, v_d_1[ k+1 ] ) ;

			y = x-k ;
			while( x < len_b && y < len_a && match( seq_b[x], seq_a[y] ) ) ++x, ++y ;
			v_d[ k ] = x ;

			if(
					(mode == myers_align_is_prefix || y == len_a) &&
					(mode == myers_align_has_prefix || x == len_b) )
			{
				char *out_a = bt_a + len_a + d +2 ;
				char *out_b = bt_b + len_b + d +2 ;
				*--out_a = 0 ;
				*--out_a = 0 ;
				for( dd = d ; dd != 0 ; )
				{
					if( k != -dd && k != dd && x == vee[ dd-1 ][ k + dd-1 ]+1 )
					{
						--dd ;
						--x ;
						--y ;
						*--out_b = seq_b[x] ;
						*--out_a = seq_a[y] ;
					}
					else if( k > -dd+1 && x == vee[ dd-1 ][ k-1 + dd-1 ]+1 )
					{
						--x ;
						--k ;
						--dd ;
						*--out_b = seq_b[x] ;
						*--out_a = '-' ;
					}
					else if( k < dd-1 && x == vee[ dd-1 ][ k+1 + dd-1 ] )
					{
						++k ;
						--y ;
						--dd ;
						*--out_b = '-' ;
						*--out_a = seq_a[y] ;
					}
					else // this better had been a match...
					{
						--x ;
						--y ;
						*--out_b = seq_b[x] ;
						*--out_a = seq_a[y] ;
					}
				}
				while( x > 0 )
				{
					--x ;
					*--out_b = seq_b[x] ;
					*--out_a = seq_a[x] ;
				}
				memmove( bt_a, out_a, bt_a + len_a + d + 2 - out_a ) ;
				memmove( bt_b, out_b, bt_b + len_b + d + 2 - out_b ) ;
				r = d ;
				goto cleanup ;
			}
		}
	}

cleanup:
	for( dd = maxd ; dd != 0 ; --dd )
		free( vee[dd-1] ) ;
	free( vee ) ;
	return r ;
}

