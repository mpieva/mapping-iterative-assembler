#ifndef INCLUDED_IMA_H
#define INCLUDED_IMA_H

#include "map_align.h"
#include "io.h"
#include "map_alignment.h"
#include "fsdb.h"
#include "pssm.h"
#include "kmer.h"
#include "assert.h"





/* init_sized_map_alignment
   Args: MapAlignmentP maln - source maln to get size from
   Returns: pointer to a fresh MapAligment for holding the
   culled results from the source maln. This culled guy
   gets the same ref and an AlnSeqArray big enough for the
   results, but no actual memory is malloced for new AlnSeq's
   Instead, we'll just point the source maln guys to this
   one if they are unique, which is determined elsewhere
*/
MapAlignmentP init_culled_map_alignment( MapAlignmentP src_maln ) ;

/* cull_maln_from_fsdb
   Args: (1) MapAlignmentP culled_maln - maln with enough room to put the
	     unique AlnSeq's
	 (2) FSDB fsdb - has valid data in front_asp, back_asp, and
	     unique_best fields
   Returns: void
   Goes through each FragSeq pointed to by fsdb->fss. For all guys that
   are unique_best and score >= SCORE_CUTOFF, copies front_asp and
   back_asp into culled_maln->AlnSeqArray. Then res
*/
void cull_maln_from_fsdb( MapAlignmentP culled_maln,
			  FSDB fsdb, int Hard_cut,
			  int SCORE_CUT_SET, double s, double n ) ;


/* asp_len
   Args: (1) AlnSeqP asp - pointer to an AlnSeq
   Returns (1) int - total length of sequence
   This function finds the total length of the sequence in this
   aligned sequence fragment. This is the sum of the sequence
   in the asp->seq field and all of the inserted sequence (if any)
   in the asp->ins array
*/
int asp_len( AlnSeqP asp ) ;


/* consensus_assembly_string
   Takes a maln object
   Generates the consensus sequence string using the aligned data
   within the maln according to the maln->cons_code and puts it
   in char* cons
   Returns char* pointer to consensus string
*/
char* consensus_assembly_string ( MapAlignmentP maln ) ;


/* Takes a pointer to an Alignment that has valid values
   in its dynamic programming matrix and valid values
   for a->aer and a->aec (ending row and column).
   Tracks back to the beginning of the alignment and
   adds valid values from a->abr and a->abc
   Returns nothing
*/
void find_align_begin( AlignmentP a ) ;



void make_ref_upper( RefSeqP ref ) ;

/* Takes a pointer to a RefSeq that has a valid
   sequence in it. Adds INIT_ALN_SEQ_LEN sequence
   from the beginning to the end so that any
   sequence fragment aligned to it will have a
   valid chance to align, despite the circularity
   of the sequence. Returns TRUE if this operation
   went smoothly, FALSE otherwise */
int add_ref_wrap( RefSeqP ref );

/* init_dpm
   Args: (1) size1 - the number of rows (fragment sequence)
         (2) size2 - the number of columns (referense sequence)
   Returns: DPMP => pointer to a dynamic programming matrix
   with memory properly allocated
*/
DPMP init_dpm( int size1, int size2 ) ;

void free_dpm( DPMP m ) ;


/* Takes a pointer to an Alignment
   that has valid sequence, length, submat, and sg data
   Does dynamic programming, filling in values in the
   a->m dynamic programming matrix
   Returns nothing */
void dyn_prog( AlignmentP a ) ;

/* size1 is length of fragment
   size2 is length of reference (wrapped if necessary) + INIT_ALN_SEQ_LEN
   rc is boolean to seay if its reverse complement
   hp_special is boolean to say if homopolymer special gap costs are to be used
*/
AlignmentP init_alignment( int size1, int size2,
			   int rc, int hp_special ) ;

/* pop_s1c_in_a
   Args: (1) AlignmentP a - has a->seq1 and a->len1 set to
   valid values
   Returns: void
   Populates the a->s1c array with code for quick lookup
   in submat
*/
void pop_s1c_in_a ( AlignmentP a ) ;

/* hp_discount_penalty
   Args: (1) int gap_len - gap length
         (2) int hplen1 - homopolymer length in column
	     (usu. longer) sequence
         (3) int hplen2 - homopolymer length in row
	     (usu. shorter) sequence
   Returns: int
   Calculates the discounted gap penalty for a gap of the input
   size that would be necessary to gap out the remaining sequence
   in the homopolymers of size indicated.
   The undiscounted (normal) gap penalty =  GOP + (GEP * length)
*/
int hp_discount_penalty ( int gap_len, int hplen1, int hplen2 ) ;

/* pop_hpl_and_hps
   Args: (1) char* seq - sequence
         (2) int len - length of sequence
	 (3) int* hpl - array of ints to be filled with hp lengths
         (4) int* hps - array of ints to be filled with hp starts
   Returns: void
   Goes through the input seq and puts in hpl the length of the
   homopolyer at each corresponding position in the sequence and
   in hps the start position of the current homopolymer. For example,
   seq=ACCGTGGTAC, len=12
   hpl=1221122111
   hps=0113455789
*/
void pop_hpl_and_hps ( const char* seq, int len, int* hpl, int* hps ) ;

/* pop_s2c_in_a
   Args: (1) AlignmentP a - has a->seq2 and a->len2 set to
   valid values
   Returns: void
   Populates the a->s2c array with code for quick lookup
   in submat
*/
void pop_s2c_in_a ( AlignmentP a ) ;


/* Input is a pointer to a valid Alignment. The value
   in a->len1 must be valid.
   Searches the last row (a->len1 -1) along all columns
   to find the best score that aligns all of the a->seq2
   Sets a->aec and a->aer to the correct values and sets
   a->best_score to the best score.
   Returns the best score */
int max_sg_score ( AlignmentP a ) ;

/* trim_frag
   Args: (1) FragSeqP frag_seq pointer to a FragSeq
         (2) char* adapter pointer to a string of the adapter
	     that may need to be trimmed
	 (3) AlignmentP align pointer to an Alignment big
	     enough for aligning the adapter to this FragSeq
   Returns: void
   Does a semi/semi global alignment and of the adapter to
   the FragSeq. If the score of this alignment is good enough
   as defined by TRIM_SCORE_CUT or if there is a perfect
   match of any number of bases >=1 at the end, then those
   are "trimmed" by setting the frag_seq->trimmed flag to
   true and the frag_seq->trim_point to the correct value
*/
void trim_frag (FragSeqP frag_seq, char* adapter,
		AlignmentP align) ;

/* Takes pointers to two PWAlnFrag's
   The first one (front_pwaln) is populated by an alignment that
   crosses the wrap_point.
   Moves all of the alignment that is behind the wrap point into
   the back_pwaln and copies over all the other info
   Sets the correct segment flag for front and back */
void split_pwaln (PWAlnFragP front_pwaln, PWAlnFragP back_pwaln,
		  int wrap_point ) ;

int populate_pwaln_to_begin( AlignmentP a, PWAlnFragP pwaln ) ;


int sg_align ( MapAlignmentP maln, FragSeqP fs, FSDB fsdb,
	       AlignmentP fw_a, AlignmentP rc_a,
	       PWAlnFragP front_pwaln,
	       PWAlnFragP back_pwaln) ;
#endif
