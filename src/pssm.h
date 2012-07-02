/* 
 * File:   pssm.h
 * Author: TCO
 *
 * Created on 3. Februar 2009, 13:14
 */

#ifndef _PSSM_H
#define	_PSSM_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "types.h"

    /* revcom_submat
       Takes a PSSMP (sm) pointer to a valid submat
       Makes a reverse complement of this submat
       Returns a pointer to the new submat
     */
    PSSMP revcom_submat(PSSMP psm);

    /* Reads in this hardcoded flat substitution matrix */
    PSSMP init_flatsubmat(void);

    /* s1b is the reference base
       s2b is the fragment (ancient) base */
    int sub_mat_score(const short int s1i,
            const short int s2i,
            int sm[][5][5],
            const int row,
            const int len);


    /* find_sm_depth
       Args: (1) int row - the current row we're on for alignment, i.e.
                           the position in the fragment sequence
             (2) len len - the length of the fragment sequence we are
                           aligning
       Returns: int - the depth in the substitution matrix for this position
                in the fragment
     */
    int find_sm_depth(int row, int len);


#ifdef	__cplusplus
}
#endif

#endif	/* _PSSM_H */

