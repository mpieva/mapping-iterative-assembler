/* 
 * File:   kmer.h
 * Author: TCO
 *
 * Created on 3. Februar 2009, 13:11
 */

#ifndef _KMER_H
#define	_KMER_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "types.h"
#include <stdlib.h>
#include <stdio.h>
#include "map_align.h"


/* add_kmer
   Args: (1) KPL* kmer array
         (2) index position - must be valid
         (3) position of this kmer to add
   Returns: void
   Takes a newly discovered kmer position and adds it to the
   array. The given index specifies what the kmer is, but
   for this operation, we really do not care. We simply add
   the position to the positions field of this kmer. If this
   kmer has never been seen before, then we have to
   initialize it, too.
*/
void add_kmer( KPL* kpa, const size_t inx, const size_t i ) ;

/* init_kpa
   Args: (1) length of kmers to use
   Returns: pointer to KPL; an array of pointers to KmerPosList
*/
KPL* init_kpa( const int kmer_len ) ;


void grow_kmers ( KmersP k ) ;

/* populate_kpa
 */
int populate_kpa( KPL* kpa, const char* seq,
		  const size_t seq_len,
		  const int kmer_len,
		  const int soft_mask ) ;

/* pop_kmers
   Args: (1) RefSeqP ref - reference sequence with forward and reverse sequence
         (2) int kmer_filt_len - length of kmers
   Initializes a Kmers struct and populates it with all the kmers in the
   forward and reverse-complement sequence of the input RefSeq.
   Returns: pointer to KmersP
*/
KmersP pop_kmers( RefSeqP ref, int kmer_filt_len ) ;

/* kmer2inx
   Args: (1) a pointer to a character string;
             the kmer to find the corresponding index of;
	     might not be null-terminated
	 (2) length of the kmer
	 (3) pointer to size_t to put the index
   Returns: TRUE if the index was set, FALSE if it could not
            be set because of some non A,C,G,T character
   Uses the formula A=>00, C=>01, G=>11, T=>11 to make a
   bit string for the kmer. Any other character is not allowed
   and will cause an error
   The bit string is constructed by reading the kmer from left
   to right. This bit-string is then interpreted as a variable
   of type size_t and is appropriate as an array index
*/
int kmer2inx( const char* kmer,
		     const unsigned int kmer_len,
		     size_t* inx ) ;

/* Returns: TRUE (1) if we should align this sequence
            FALSE (0) if we should NOT align this sequence because
                      it shares no kmers with the reference
*/
int new_kmer_filter( FragSeqP fs,
		     KPL* fkpa,
		     KPL* rkpa,
		     int kmer_len,
		     AlignmentP fwa,
		     AlignmentP rca ) ;

int kmer_filter( int kmer_filt_len, FragSeqP fs, KmersP k ) ;


#ifdef	__cplusplus
}
#endif

#endif	/* _KMER_H */

