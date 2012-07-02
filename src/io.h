/* 
 * File:   io.h
 * Author: Ed Green
 *         Michael Siebauer
 *
 * Created on 25. Januar 2009, 14:41
 */

#ifndef _IO_H
#define	_IO_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "params.h"
#include "types.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "map_align.h"

/* find_input_type
   Args: 1. FILE* pointer to file to be analyzed
   Returns: sequence code indicating what kind of sequence file
            this is:
	    0 => fasta
	    1 => fastq
   Resets the input FILE pointer to the beginning of the file
*/
  int find_input_type( FILE * FF );

/* read_next_seq
   Args: 1. FILE* pointer to file being read
         2. FragSeqP pointer to FragSeq where the next sequence data will go
	 3. int code indicating which parser to use
   Returns: TRUE if a sequence was read,
            FALSE if EOF
*/

  int read_next_seq( FILE * FF, FragSeqP frag_seq, int seq_code );

/* read_fasta
   args 1. pointer to file to be read
        2. pointer to FragSeq to put the sequence
   returns: TRUE if sequence was read,
            FALSE if EOF or not fasta
*/
int read_fasta ( FILE * fasta, FragSeqP frag_seq );

/* read_fastq
   Args 1. pointer to file to be read
        2. pointer to FragSeq to put the sequence into
   Returns: TRUE if a sequence was read,
            FALSE if EOF
*/

int read_fastq ( FILE * fastq, FragSeqP frag_seq );

/* calc_qual_sum
   Args: 1. pointer to a string of quality scores for this sequence
   Returns: 1. int - the sum of quality scores for this sequence
   This assumes that quality scores are represented as the 
   ASCII code + 64
*/
  inline int calc_qual_sum( const char* qual_str );



/* Read in the reference sequence from fasta file and make reverse complement, too
 * Return 1 success
 0 failure
*/
int read_fasta_ref(RefSeqP ref, const char* fn);


  /* Reads in a set of scoring matrices for each of the PSSM_DEPTH
     positions at the beginning and end of the sequence that are
     to have special scoring matrices and the single 'MIDDLE' matrix
     for everything in the middle.
     Puts these matrices into a PSSM structure.
     Returns a pointer to this structure (PSSMP)
  */
  PSSMP read_pssm(const char* fn);
  
  /* Reads one pairwise alignment from an Udo Stenzel align
     output file of semi-global alignments against a common
     target sequence (usually chrM) into a PWAlnFrag.
     Args: FILE* advanced to next pairwise alignment
     PWAlnFragP to be populated
     Returns 1 if success;
     0 if EOF or failure
     -1 for failure
     */
  int read_align_aln(FILE* align_f, PWAlnFragP af);
  
  FILE * fileOpen(const char *name, char access_mode[]);
  
  IDsListP parse_ids(char* fn);
  
  // Prints this string colored
void color_print(char* string);


  void ace_output(MapAlignmentP maln);

  void line_print_cons(char* consensus, char* aln_ref, char* ref_id, int* cov);

  void clustalw_print_cons(char* cons, char* aln_ref, char* ref_id);



#ifdef	__cplusplus
}
#endif

#endif	/* _IO_H */

