
#ifndef INCLUDED_MAP_ALIGN_H
#define INCLUDED_MAP_ALIGN_H

#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <getopt.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include "params.h"
#include "types.h"
#include <string.h>
#include "io.h"
#include "map_alignment.h"


/* Function Prototypes */


/* base2inx
 Args (1) char base - the base to find the index for
 Returns a short int of the index position into a substitution
 matrix for this base
 This function finds the corresponding index for a
 PSSMP->sm[DEPTH][][] (position-specific substitution matrix)
 for a given base. This index is used to look up the appropriate
 score. This function assumes that the rows and columns are in
 the order: A, C, G, T, N
 */

inline short int base2inx(const char base) ;

int idCmp(const void* id1_, const void* id2_) ;

/* Takes an AlnSeqP and a beginning and end coordinate of a region.
 All coordinates are 0-based
 Returns true is this AlnSeq overlaps the region at all, false
 if it does not */
inline int alnseq_ol_reg(AlnSeqP as, const int rs, const int re);


/* This IDsList */
IDsListP init_ids_list(void) ;

void add_id(char* new_id, IDsListP used_ids_list) ;

void grow_ids_list(IDsListP ids) ;

int allowed_alignment(int ids_rest, IDsListP rest_ids_list, int no_dups,
		IDsListP used_ids_list, PWAlnFragP pwaln, double score_int,
		double score_slo) ;

void show_single_pos(int ref_pos, char ref_base, char cons_base, BaseCountsP bcs) ;

void add_base(char b, BaseCountsP bcs, PSSMP psm, int pssm_code) ;

void reset_base_counts(BaseCountsP bc) ;

/* Takes a pointer to a BaseCounts bcs and the maln->cons_code
   (consensus code)
   The bcs must have valide data
   Returns the character of the consensus at this position as
   defined by the scheme to use and the data in bcs
   If the coverage is 0, returns N
   If there are >= PERC4GAP percent of reads with a gap, returns -
   Otherwise:
   If cons_code = 1, returns any base with score >= MIN_SCORE_CONS
                     or N if none
      cons_code = 2, returns any base with score >= MIN_SC_DIFF_CONS
                     better than the second highest scoring base or
		     N if none
*/
char find_consensus(BaseCountsP bcs, int cons_code) ;

int alnSeqCmp(const void* as1_, const void* as2_) ;

char revcom_char(const char base) ;

/* Takes a MapAlignmentP and a position where some of
 the aligned fragments have an insert relative to the
 reference. That is, maln->ref->gaps[position] > 0.
 Populates the char* ins_cons and int* cons_cov
 arrays with the consensus sequence and consensus
 coverage, respectively. These must be appropriately
 sized elsewhere. If out_format is the special value
 of 4, then we just show these differences now and
 do not return anything.
 */
void find_ins_cons(MapAlignmentP maln, int pos, char* ins_cons, int* cons_cov,
		int out_format) ;

void revcom_PWAF(PWAlnFragP pwaln) ;

/* For a given region, defined by reg_start and reg_end, show
 the refence sequence, the consensus sequence,
 and the sequence of all the fragments that overlap this
 region at all.
 */
void print_region(MapAlignmentP maln, int reg_start, int reg_end,
		int out_format, int in_color) ;

void col_print_cons(char* consensus, char* aln_ref, int* cov, int* ref_poss,
		MapAlignmentP maln) ;

/* Takes a pointer to a populated PWAlnFrag (pwaln) and
 a pointer to a populated MapAlignent (maln)
 Does:
 1. Adds this aligned sequence, without gaps to maln->AlnSeqArray,
 growing this array if necessary
 2. Populates the gaps array of this newly aligned fragment to
 indicate where its gaps are relative to the reference
 3. Updates the gaps array of the reference sequence (maln->ref->gaps[])
 and the gaps array of all aligned fragments to accomodate any
 new gaps this new fragment may require
 Returns: 1 (TRUE) if success
 0 (FALSE) if failure
 */
int merge_pwaln_into_maln(PWAlnFragP pwaln, MapAlignmentP maln) ;


/* Takes the description from an Udo align aligned sequence
 and puts the start, end, strand, and score information in
 the correct field of the PWAlnFragP
 The desc (description) is a string like this, e.g.:
 "- 4199-4261 score=5441"
 */
int ses_from_align_desc(PWAlnFragP pwaln, int* strand) ;

/* adapt_from_desc checks the frag_desc string of a PWAlnFrag,
 given a pointer to one (PWAlnFragP) and sets the trimmed
 field to true (1) if the phrase "adapter cut off" is there
 Returns true if everthing went fine, false otherwise
 */
int adapt_from_desc(PWAlnFragP af) ;

/* Grow the space for a sequence (an array of char)
 to twice its current size
 Copy its current contents into the new sequence
 Free the now unused old memory
 */
char* grow_seq(char* seq, int size);





#endif
