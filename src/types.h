/*
 * File:   types.h
 * Author: 
 *
 * Created on 16. Januar 2009, 14:59
 */

#ifndef _TYPES_H
#define	_TYPES_H

#include "params.h"
#include <stdlib.h>
#include <ctype.h>


#define save_malloc malloc

#ifdef	__cplusplus
extern "C" {
#endif




/*
  Define PWAlnFrag as a struct pw_aln_frag
  This is a structure for holding an alignment pair
  One of the pair is an aligned fragment of the
  reference sequence. The other is a sequence fragment
  aligned to it.
  This is used as a temporary holding place for alignments
  until they're merged into the big MapAlignment
  There is no functionality to grow this guy, so it should
  be as big as necessary from the beginning. This is
  determined by INIT_ALN_SEQ_LEN
*/
typedef struct pw_aln_frag {
  char ref_id[MAX_ID_LEN + 1];
  char ref_desc[MAX_DESC_LEN + 1];
  char frag_id[MAX_ID_LEN + 1];
  char frag_desc[MAX_DESC_LEN + 1];
  char ref_seq[(2*INIT_ALN_SEQ_LEN)+1];
  char frag_seq[(2*INIT_ALN_SEQ_LEN)+1];
  int start;
  int end;
  int revcom;
  int trimmed;
  int score;
  char segment; // f=front, a=all, b=back, n=not applicable
  int num_inputs; // for collapsed sequences, the number of input seqs
  int offset; // for segment='b', number of bases not shown that
              // were in the front fragment
} PWAlnFrag;
typedef struct pw_aln_frag* PWAlnFragP;

/*
   Define Alnseq as a struct aln_seq
   This is simply a structure for holding a
   string of aligned sequence
 */
typedef struct alnseq {
  char id[MAX_ID_LEN + 1]; // the ID of the sequence
  char desc[MAX_DESC_LEN + 1]; // the description of the sequence
  char seq[ (2*INIT_ALN_SEQ_LEN) + 1];  // the sequence string
  char smp[ (2*INIT_ALN_SEQ_LEN) + 1];  // code for substitution matrix depth
  char* ins[ (2*INIT_ALN_SEQ_LEN) + 1]; // array of pointers to char
  // that will be filled with sequence
  int start;  // where this sequence starts relative to the reference (0-indexed)
  int end;    // where this sequence ends relative to the reference (0-indexed)
  int revcom; // boolean to denote that this sequence has been
              // reverse complemented
  int trimmed; // boolean to denote that this sequence has been trimmed
  int score;  // the alignment score for this guy
  int num_inputs; // the number of input seqs if this is a collapsed seq
  char segment; // f=front, a=all, b=back, n=not applicable
} AlnSeq;
// pointer to struct aln_seq
typedef struct alnseq* AlnSeqP;

/*
  Define RefSeq and RefSeqP to be the reference sequence
  against which all the fragments have been aligned.
*/
typedef struct refseq {
  char id[MAX_ID_LEN + 1]; // the ID of the reference sequence
  char desc[MAX_DESC_LEN + 1]; // the description of the sequence
  char* seq;                   // the sequence as a string
  char* rcseq;                 // the reverse complement as a string
  int seq_len;                 // the length of the sequence
  int size;                    // size of char array for this aligned sequence
  int* gaps;               // array giving the size of the longest gap
                           // that has been introduced at each position
  int circular;            // Boolean to denote circular sequence
  int wrap_seq_len;        // length of sequence with extra wrapped bit
  // seq_len remains the actual length of the sequence
} RefSeq;
// pointer to struct refseq
typedef struct refseq* RefSeqP;

/* Define qsumseq */
typedef struct qsumseq {
  unsigned int Aqualsum[INIT_ALN_SEQ_LEN+1];
  unsigned int Cqualsum[INIT_ALN_SEQ_LEN+1];
  unsigned int Gqualsum[INIT_ALN_SEQ_LEN+1];
  unsigned int Tqualsum[INIT_ALN_SEQ_LEN+1];
} QSumSeq;
typedef struct qsumseq* QSSP;

/* Define FragSeq and FragSeqP to hold a simple sequence */
typedef struct fragseq {
  char id[MAX_ID_LEN + 1];
  char desc[MAX_DESC_LEN + 1];
  char seq[INIT_ALN_SEQ_LEN+1];
  char qual[INIT_ALN_SEQ_LEN+1];
  QSSP qss; // pointer to a QSumSeq struct that may be needed for collapsing
  int qual_sum;
  int trim_point; // 0-indexed position of last base before adapter
  int trimmed; // Boolean, TRUE means sequence should be trimmed to trim_point
  int seq_len;
  int strand_known;  // Boolean, TRUE means the alignment strand of this
  // sequence has been learned by virtue of a positive scoring alignment
  int rc; // Boolean, TRUE means this is the reverse complement
  int as; // 0-indexed start point of alignment on current ref
  int ae; // 0-indexed end point of alignment on current ref
  int score; // current score of alignment on reference
  AlnSeqP front_asp; // pointer to where I can find the front AlnSeq
  AlnSeqP back_asp;  // pointer to where I can find the back AlnSeqP
  //                   (if applicable, otherwise NULL)
  int unique_best;   // boolean; TRUE means unique & best score
  //                    for repeat filtering
  int num_inputs; // number of sequences collapsed into this one
} FragSeq;
typedef struct fragseq* FragSeqP;

/* Define fragseqdb and FSDB to hold a database of FragSeqs */
typedef struct fragseqdb {
  FragSeqP* fss; // Pointer to array of FragSeqs
  int       trim_sort; // Use extra information about whether sequence is trimmed when
                       // determining whether a sequence is unique
  size_t    size; // Current size of array pointed to by fss
  size_t    num_fss; // Current number of FragSeqs in fss
} FragSeqDB;
typedef struct fragseqdb* FSDB;

/* Define PSSM as an array of position specific substitution
   matrices to be used at different points in an alignment.
   The first PSSM_DEPTH-1 matrices are for the beginning of the
   alignment. The matrix at PSSM_DEPTH is for the middle. The
   matrices at PSSM_DEPTH+1..2*PSSM_DEPTH are for the end of
   the alignment.
   The second value is for the reference base. The third value
   is for the ancient base.
   A=0, C=1, G=2, T=3, anthing else=4
*/
typedef struct pssm {
  int sm[2*PSSM_DEPTH+1][5][5];
  int depth;
} PSSM;
typedef struct pssm* PSSMP;

/* Define DPE to be an element of a dynamic programming
   matrix. Each element remembers its best score and
   keeps a pointer to where it came from */
typedef struct dpe {
  int score; // best score at this position
  int trace; // code for previous aligned position that
  // gave the best score; 0 => diagonal;
  // pos. number => column number
  // neg. number => row number
  // If number == current row, this is the beginning
  // of the alignment
} DPE;
typedef struct dpe* DPEP;

typedef struct dpm {
  DPEP* mat;
  int rows;
  int cols;
} Mat;
typedef struct dpm* DPMP;


typedef struct map_alignment {
  RefSeqP ref;       // The reference sequence to which everything is mapped
  PSSMP fpsm;        // The PSSMP set of + strand matrices for aligning and consensus
  PSSMP rpsm;        // The PSSMP set of - strand matrices for aligning and consensus
  int num_aln_seqs;  // Number of sequences in this alignment
  int size;          // Length of AlnSeqArray
  int cons_code;     // Code for scheme for determining the consensus base
                     //    1 => only majority rule consensus
                     //    2 => (unique) plurality rule consensus
  int distant_ref;   // initial reference sequence is distantly related
  AlnSeqP* AlnSeqArray;
} MapAlignment;
// pointer to struct alignment
typedef struct map_alignment* MapAlignmentP;

typedef struct base_counts {
  int As;
  int scoreA;
  int Cs;
  int scoreC;
  int Gs;
  int scoreG;
  int Ts;
  int scoreT;
  int gaps;
  int cov;
  double frac_agree; // fraction of bases that agree with the most
                     // common base
} BaseCounts;
typedef struct base_counts* BaseCountsP;

typedef struct alignment {
  const char* seq1; // reference sequence
  const char* seq2; // fragment sequence
  short int* s1c; // array of submat lookup indeces for s1 - must be
            // dynamically allocated
  short int s2c[INIT_ALN_SEQ_LEN]; // code for submat lookup for
                             // sequence2 that cannot be longer
  int len1;   // length of reference sequence
  int len2;   // length of fragment sequence
  unsigned char* align_mask; // 0 => alignment cannot be here;
                             // 1 => alignment can go through here

  PSSMP submat;  // position substitution matrices
  int gop;    // gap open penalty
  int gep;    // gap extension penalty
  int hp;     // Boolean, TRUE = special discount for homopolymer
              // associated gaps
  int* hpcl;  // array of lengths of hps for each seq1 position
  int* hpcs;  // array of starts of hps for each seq1 position
  int* hprl;  // array of lenghs of hps for each seq2 position
  int* hprs;  // array of starts of hps for each seq2 position
  DPMP m;     // pointer to struct dpm, dynamic prog. matrix
  int* best_gap_row; // array of current best row to gap to,
  //                     useful during dynaminc programming
  int best_gap_col; // keeps column number of current best-
  //                   scoring gap column
  int sg5;    // Boolean, TRUE = do semiglobal alignment at 5' end of
  //                             seq2 (pay penalty for unaligned)
  //                      FALSE = local alignment at 5' end
  int sg3;    // Boolean, TRUE = do semiglobal alignment at 3' end of
  //                             seq2 (pay penalty for unaligned)
  //                      FALSE = local alignment at 3' end
  int rc;     // Boolean, TRUE => this is an alignment to the
  // reverse complement, FALSE => forward strand alignment
  int abc; // alignment beginning column
  int abr; // alignment beginning row
  int aec; // alignment ending column
  int aer; // alignment ending row
  int best_score; // score at m->[aer][aec], i.e., the best score
} Alignment;
typedef struct alignment* AlignmentP;

typedef struct ids_list {
  int num_ids;
  int sorted;
  int size;
  char** ids;
} IDsList;
typedef struct ids_list* IDsListP;

typedef struct kmers {
  int num_kmers;
  int kmer_len;
  int sorted;
  int size;
  char** kmers;
} Kmers;
typedef struct kmers* KmersP;

typedef struct kmer_pos_list {
  size_t num_pos; // current number of known positions for this kmer
  int    sorted;  // boolean; TRUE => positions are in order
  unsigned int positions[MAX_KMER_POS]; // list of known positions
                                        // for this kmer
} KmerPosList;
typedef struct kmer_pos_list* KPL;





#ifdef	__cplusplus
}
#endif

#endif	/* _TYPES_H */

