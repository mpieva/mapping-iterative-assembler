/* 
 * File:   params.h
 * Author:
 *
 * Created on 16. Januar 2009, 15:02
 */

#ifndef _PARAMS_H
#define	_PARAMS_H

#ifdef	__cplusplus
extern "C" {
#endif

#define DEBUG (0)
#define CONS_SCHEME (1)
#define MAX_ID_LEN (64)
#define MAX_DESC_LEN (128)
#define CLUSTALW_LINE_WIDTH (60)
#define FASTA_LINE_WIDTH (60)
#define MAX_LINE_LEN (1000000)
#define PSSM_DEPTH (15)
#define MAX_FN_LEN (1023)
#define SCORE_CUTOFF_BUFFER (80) // just a guess for now
#define FIRST_ROUND_SCORE_CUTOFF (2000) // reference alignment original cutoff
#define GOP (1000) // Gap open penalty
#define GEP (200) // Gap extension penalty
#define FLAT_MATCH (200) // score in the flat matrix for a gap
#define FLAT_MISMATCH (-600) // score in the flat matrix for a mismatch
#define N_SCORE (-100)
#define NR_SCORE (-10) // score for N in reference
#define TRIM_SCORE_CUT (1000)
#define MAX_ITER (30) // maximum number of assembly iterations to do
#define REALIGN_BUFFER (50) // amount of sequence padding to add in realignment
#define QUAL_ASCII_OFFSET (33) // ascii code of lowest quality score, i.e. 0
#define DEF_S 200.0
#define DEF_N 0.0
#define MIN_ALIGNABLE_LEN (15) // when distant reference is used, minimum amount of
  // alignable sequence when reducing the sequence length for bases that overlap
  // N positions in the reference
#define MIN_SCORE_CONS (-399) // minimum score to call consensus base, not N, under
                              // cons_code 1
#define MIN_SC_DIFF_CONS (2400) // minimum diff between best and 2nd best to call
                                // best base consensus under cons_code 2
#define PERC4GAP 50 // minimum percent of reads at a position that have a
                    // gap for the consensus to be a gap

/* INIT_NUM_IDS is the initial number of sequence IDs that
   can be in the file that we're restricting analysis to
*/
#define INIT_NUM_IDS (1048576)

/* MAX_INS_LEN is the size of the char array accomodating
   sequence inserts in an aligned fragment relative to the
   reference sequence. That is, it's the longest single
   gaps size allowable in the reference sequence alignment
*/
#define MAX_INS_LEN (512)

/* INIT_REF_SEQ_LEN is initial size of refe::rence sequence to
   which all aligned fragments are mapped. It can grow when
   read in, if necessary */
#define INIT_REF_SEQ_LEN (32768)

/* INIT_ALN_SEQ_LEN is the initial and maxmial length of
   aligned sequence fragments. It cannot grow, so make sure
   this is big enough */
#define INIT_ALN_SEQ_LEN (256)
#define INIT_NUM_ALN_SEQS (16000)


#define MAX_FN_LEN (1023)


#define MAX_KMER_POS (128)
#define MAX_KMER_LEN (14)
#define KMER_SATURATE (128)
#define ALIGN_MASK_BUFFER (10)




#ifdef	__cplusplus
}
#endif

#endif	/* _PARAMS_H */

