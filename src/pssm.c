
#include "pssm.h"

/* s1b is the reference base
   s2b is the fragment (ancient) base */
int sub_mat_score(const short int s1i,
        const short int s2i,
        int sm[][5][5],
        const int row,
        const int len) {
    int score;

    if (row < PSSM_DEPTH) {
        score = sm[row][s1i][s2i];
        return score;
    }

    if (len - (row + 1) < PSSM_DEPTH) {
        score = sm[(PSSM_DEPTH * 2)-(len - (row + 1))][s1i][s2i];
        return score;
    }

    score = sm[PSSM_DEPTH][s1i][s2i];
    return score;

}

/* find_sm_depth
   Args: (1) int row - the current row we're on for alignment, i.e.
                       the position in the fragment sequence
         (2) len len - the length of the fragment sequence we are
                       aligning
   Returns: int - the depth in the substitution matrix for this position
            in the fragment
 */
int find_sm_depth(int row, int len) {
    if (row < PSSM_DEPTH) {
        return row;
    }

    if (len - (row + 1) < PSSM_DEPTH) {
        return ( (PSSM_DEPTH * 2)-(len - (row + 1)));
    }

    return PSSM_DEPTH;
}

/* revcom_submat
   Takes a PSSMP (sm) pointer to a valid submat
   Makes a reverse complement of this submat
   Returns a pointer to the new submat
 */
PSSMP revcom_submat(PSSMP psm) {
    int d, rcd;
    PSSMP rcpsm;
    rcpsm = (PSSMP) save_malloc(sizeof (PSSM));
    rcpsm->depth = PSSM_DEPTH;

    for (d = 0; d <= (PSSM_DEPTH * 2); d++) {
        rcd = (PSSM_DEPTH * 2) - d;

        rcpsm->sm[rcd][0][0] = psm->sm[d][3][3];
        rcpsm->sm[rcd][0][1] = psm->sm[d][3][2];
        rcpsm->sm[rcd][0][2] = psm->sm[d][3][1];
        rcpsm->sm[rcd][0][3] = psm->sm[d][3][0];
        rcpsm->sm[rcd][0][4] = psm->sm[d][3][4];

        rcpsm->sm[rcd][1][0] = psm->sm[d][2][3];
        rcpsm->sm[rcd][1][1] = psm->sm[d][2][2];
        rcpsm->sm[rcd][1][2] = psm->sm[d][2][1];
        rcpsm->sm[rcd][1][3] = psm->sm[d][2][0];
        rcpsm->sm[rcd][1][4] = psm->sm[d][2][4];

        rcpsm->sm[rcd][2][0] = psm->sm[d][1][3];
        rcpsm->sm[rcd][2][1] = psm->sm[d][1][2];
        rcpsm->sm[rcd][2][2] = psm->sm[d][1][1];
        rcpsm->sm[rcd][2][3] = psm->sm[d][1][0];
        rcpsm->sm[rcd][2][4] = psm->sm[d][1][4];

        rcpsm->sm[rcd][3][0] = psm->sm[d][0][3];
        rcpsm->sm[rcd][3][1] = psm->sm[d][0][2];
        rcpsm->sm[rcd][3][2] = psm->sm[d][0][1];
        rcpsm->sm[rcd][3][3] = psm->sm[d][0][0];
        rcpsm->sm[rcd][3][4] = psm->sm[d][0][4];

        rcpsm->sm[rcd][4][0] = psm->sm[d][4][3];
        rcpsm->sm[rcd][4][1] = psm->sm[d][4][2];
        rcpsm->sm[rcd][4][2] = psm->sm[d][4][1];
        rcpsm->sm[rcd][4][3] = psm->sm[d][4][0];
        rcpsm->sm[rcd][4][4] = psm->sm[d][4][4];
    }
    return rcpsm;
}

/* Reads in this hardcoded flat substitution matrix */
PSSMP init_flatsubmat(void) {
    PSSMP flatsubmat;
    int cur_pos, base, other_base;
    flatsubmat = (PSSMP) save_malloc(sizeof (PSSM));
    flatsubmat->depth = PSSM_DEPTH;

    for (cur_pos = 0; cur_pos <= (2 * PSSM_DEPTH); cur_pos++) {
      for (base = 0; base <= 4; base++) {
	for (other_base = 0; other_base <= 3; other_base++) {
	  if (base == other_base) {
	    flatsubmat->sm[cur_pos][base][other_base]
	      = FLAT_MATCH;
	  } 
	  else {
	    flatsubmat->sm[cur_pos][base][other_base]
	      = FLAT_MISMATCH;
	  }
	}
	// last column is for N
	flatsubmat->sm[cur_pos][base][4] = N_SCORE;
      }
    }

    /* Special score, NR_SCORE, for reference having an N */
    for ( cur_pos = 0; cur_pos <= (2*PSSM_DEPTH); cur_pos++ ) {
      for ( base = 0; base <= 4; base++ ) {
	flatsubmat->sm[cur_pos][4][base] = NR_SCORE;
      }
    }
    return flatsubmat;
}


