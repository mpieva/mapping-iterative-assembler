/* 
 * File:   fsdb.h
 * Author: TCO
 *
 * Created on 3. Februar 2009, 13:08
 */

#ifndef _FSDB_H
#define	_FSDB_H

#include "types.h"
#include "params.h"
#include "stdio.h"
#include "stdlib.h"
#include "io.h"

#ifdef	__cplusplus
extern "C" {
#endif



/* fs_comp
   Args: (1) pointer to first FragSeqP
         (2) pointer to second FragSeqP
   Returns: 1 if the first FragSeq comes before the second one,
           -1 if it comes after and
	   0 if they come at the same time
   This function defined a sort order for FragSeqP's. This order
   is useful for then determining which FragSeqP's are unique.
*/
  int fs_comp ( const void* fs1_,
		const void* fs2_ ) ;

  int fs_comp_qscore ( const void* fs1_,
		       const void* fs2_ );
/* add_virgin_fs2fsdb
   Args: (1) FragSeqP fs - pointer to a "virgin" FragSeq
         (2) FSDB fsdb - database to add this FragSeq to
   Returns: 1 if success; 0 if failue (not enough memories)
   This function is only called from sg_align; the argument
   FragSeqP points to a FragSeq for which the following is
   true: id, desc, as, ae, score, front_asp, back_asp and
   unique are set to correct values.
   If trimmed is true, then this sequence is to be trimmed
   to the trim_point
   If rc is set, then this sequence is to be reverse
   complemented
   Once these operations are done, this "non-virgin" FragSeq
   is then copied into the next slow of fsdb, growing fsdb
   if necessary, and incrementing its fsdb->num_fss
*/
  int add_virgin_fs2fsdb( FragSeqP fs, FSDB fsdb ) ;

/* Sorts the fsdb->fss on rc, as, ae, score
   After sorting all 1 strand alignments are first
   These are sorted by as, then ae, with the highest
   scoring guys first and then lower scoring guys
*/
  void sort_fsdb( FSDB fsdb ) ;
  void sort_fsdb_qscore( FSDB fsdb );


/* find_fsdb_score_cut
   Args: (1) FSDB fsdb - has valid data for seq_len, score, and
             unique_best
     (2) double* slope - pointer to slope to be calculated
     (3) double* intercept - pointer to intercept to be calc.
   Returns: void
   Takes all the sequence lengths and scores in a FSDB database.
   Calculates the best fit line through the data:
   score = (slope * seq_len) + intercept
   That is, is determines the dependency of average score on the
   length of the sequence. This can then be used to determine
   what is an inappropriately scoring (for its length) alignment.
*/
  void find_fsdb_score_cut( FSDB fsdb, double* slope, double* intercept ) ;
  
/* write_fastq
   Args: (1) char* fn
         (2) FSDB fsdb
   Returns: void
   Writes a fastq database of sequences to the filename given of all
   sequences and quality scores in the fsdb
*/
  void write_fastq( char* fn, FSDB fsdb );


/* set_uniq_in_fsdb
   Args: (1) FSDB fsdb - has fss field SORTED!
         (2) int just_outer_coords - boolean; TRUE means just use
	     outer coordinate info (strand, start, end) to
	     decide about uniqueness; FALSE is a more complex
	     scheme. If a sequence is has the same start, but a
	     lower end point, it is not unique unless it is also
	     trimmed. This is to handle 454 data where occassionally
	     sequences "end" because of some filter, but not the
	     natural end of the molecule. Then, repeats can show up
	     in different lengths!
   Returns: void
   Goes through each sequence and the sets the unique_best flag
   to true for the first of each kind (same as, ae, and rc) and
   sets unique_best to false for all others
*/

  void set_uniq_in_fsdb( FSDB fsdb, const int just_outer_coords ) ;

/* pop_smp_from_FSDB
   Args: (1) FSDB fsdb with valid data
         (2) Depth of PSSM matrices
   Returns: void
   Goes through all the AlnSeqs in the fsdb->fss array. Follows
   the front_asp (and back_asp, if necessary) pointer to populate
   the smp field of all AlnSeqs with the correct code for what
   depth in the PSSM matrix to use for consensus calling */
void pop_smp_from_FSDB( FSDB fsdb, int depth ) ;

/* add_fs2fsdb
   Args: (1) FragSeqP fs - pointer to a fully valid FragSeq
         (2) FSdb fsdb - database to add this FragSeq to
   Returns: 1 if success; 0 if failure (not enough memories)
   Adds the FragSeq pointed to by fs to the fsdb database,
   growing it if necessary */
int add_fs2fsdb( FragSeqP fs, FSDB fsdb ) ;

/* grow_FSDB
   Args: (1) FSDB (fsdb) to be made twice as big
   Returns: 1 if success; 0 if failure (not enough memories)
   Grows an FSDB by allocating another chunk of memory for
   the FragSeqs as big as the one it already has. Note, it
   *DOES NOT* throw away the one it already has. Then, the
   fsdb->fss array is replaced by one twice as big. The
   pointers to all the existing FragSeqs are copied over
   and the new ones are set up. The size is reset, too.
   The old fsdb->fss array is freed
*/
int grow_FSDB( FSDB fsdb ) ;

/* init_FSDB
   Arguments: void
   Returns: FSDB (pointer to struct fragseqdb) / NULL if
    not enough memories
   Used for initializing a new database of FragSeqs. Allocates
   enough memoery for INIT_NUM_ALN_SEQS of these
*/
FSDB init_FSDB ( void );




#ifdef	__cplusplus
}
#endif

#endif	/* _FSDB_H */

