/*
 * File:   map_alignment.h
 * Author: TCO
 *
 * Created on 26. Januar 2009, 12:16
 */

#ifndef _MAP_ALIGNMENT_H
#define	_MAP_ALIGNMENT_H

#include "types.h"
#include "io.h"
#include "config.h"

#ifdef	__cplusplus
extern "C" {
#endif



    /* Initialize a MapAlignment object and return a pointer to it */
    MapAlignmentP init_map_alignment(void);



    /* free_map_alignment
     Takes a MapAlignmentP (maln)
     Frees the memory pointed to by its components
     Returns nothing
     */
    void free_map_alignment(MapAlignmentP maln);


    /* Write out the data in a MapAlignment data structure
     to a file
     */
    int write_ma(char* fn, MapAlignmentP maln);

    MapAlignmentP read_ma(const char* fn);

    /* Grow the space for a MapAlignment to twice its current
 size. Actually, just grow the array of aligned sequences.
 Copy the current aligned sequences into the new array
 free the now unused memory.
 Return 1 if success
 0 if failure
     */
    int grow_alns_map_alignment(MapAlignmentP aln);

    // Return the absolute number of gaps upstream of this position
    int sum_of_gaps(MapAlignmentP maln, int pos);


    MapAlignmentP init_map_alignment(void);
    int count_aln_seqs(MapAlignmentP maln);
    void sort_aln_frags(MapAlignmentP maln);




    void show_consensus(MapAlignmentP maln, int out_format);




    //some more convenient functions
    int get_consensus_length(MapAlignmentP maln);
    char* get_consensus(MapAlignmentP maln);


    void print_assembly_summary(MapAlignmentP maln);

    /* Sorts the AlnSeqArray by alnSeqCmp (its start and end
     coordinates.
     Note that after this operation, any FragSeqDB pointing
     to this AlnSeqArray will be wrong! */
    void sort_aln_frags(MapAlignmentP maln);


#ifdef	__cplusplus
}
#endif

#endif	/* _MAP_ALIGNMENT_H */

