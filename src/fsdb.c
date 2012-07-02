#include "fsdb.h"



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
	      const void* fs2_ ) {
  FragSeqP* fs1p = (FragSeqP*) fs1_;
  FragSeqP* fs2p = (FragSeqP*) fs2_;
  FragSeqP fs1 = *fs1p;
  FragSeqP fs2 = *fs2p;

  /* First sort criteria is strand (rc) */
  if ( fs1->rc && !(fs2->rc) ) {
    return -1;
  }
  if ( !(fs1->rc) && fs2->rc ) {
    return 1;
  }

  /* Forward strand guys */
  if ( fs1->rc == 0 ) {
    /* Second sort criteria is where they start,
       lower coordinates come first */
    if ( fs1->as < fs2->as ) {
      return -1;
    }
    if ( fs1->as > fs2->as ) {
      return 1;
    }
    /* Third sort criteria is where they end,
       lower coordinates come later */
    if ( fs1->ae < fs2->ae ) {
      return 1;
    }
    if ( fs1->ae > fs2->ae ) {
      return -1;
    }
    /* Fourth sort criteria is the score,
       lower score comes later */
    if ( fs1->score < fs2->score ) {
      return 1;
    }
    if ( fs1->score > fs2->score ) {
      return -1;
    }

    /* If they match on all that, they are the same */
    return 0;
  }

  /* Reverse strand guys */
  else {
    /* Sort by end (start of molecule, higher coordinates first */
    if ( fs1->ae < fs2->ae ) {
      return 1;
    }
    if ( fs1->ae > fs2->ae ) {
      return -1;
    }

    /* Sort by start (end of molecule, higher coordinates later */
    if ( fs1->as < fs2->as ) {
      return -1;
    }
    if ( fs1->as > fs2->as ) {
      return 1;
    }

    /* Sort by score, lower scores come later */
    if ( fs1->score < fs2->score ) {
      return 1;
    }
    if ( fs1->score > fs2->score ) {
      return -1;
    }

    /* If all that matches, they are sorted the same */
    return 0;
  }
}



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
int add_virgin_fs2fsdb( FragSeqP fs, FSDB fsdb ) {
  int i, len, half_len;
  char tmp_b;

  /* Trim it? */
  if ( fs->trimmed ) {
    fs->seq[fs->trim_point + 1] = '\0';
    fs->seq_len = fs->trim_point + 1;
  }

  /* revcom it? */
  if ( fs->rc ) {
    len = fs->seq_len;
    half_len = len / 2;
    for ( i = 0; i < half_len; i++ ) {
      tmp_b = fs->seq[i];
      fs->seq[i] = revcom_char(fs->seq[len-(i+1)]);
      fs->seq[len-(i+1)] = revcom_char(tmp_b);
    }
    if ( len%2 == 1 ) {
      /* Sequence length was odd, revcom the middle base */
      fs->seq[half_len] = revcom_char(fs->seq[half_len]);
    }
  }

  /* OK, now copy it over to fsdb */
  return ( add_fs2fsdb( fs, fsdb ) );
}



/* Sorts the fsdb->fss on rc, as, ae, score
   After sorting all 1 strand alignments are first
   These are sorted by as, then ae, with the highest
   scoring guys first
*/
void sort_fsdb( FSDB fsdb ) {
  qsort( (void*)fsdb->fss, (size_t)fsdb->num_fss,
	 sizeof(FragSeqP), fs_comp );
}


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
void find_fsdb_score_cut( FSDB fsdb, double* slope, double* intercept ) {
  /* Load up the lengths and scores of all unique guys that
     had sensible scores as defined by FIRST_ROUND_SCORE_CUTOFF
     This is necessary in case the distant reference option is
     used in which case we may have some total crap sequences
     and scores that will screw up the fit.  (This exact same loop is
     repeated twice (three times in debug mode), the single pass
     algorithm would be numerically unstable.)
  */
  size_t j = 0, i ;

  double xbar = 0, ybar = 0 ;
  for ( i = 0; i < fsdb->num_fss; i++ ) {
    if ( fsdb->fss[i]->unique_best &&
	(fsdb->fss[i]->score >= FIRST_ROUND_SCORE_CUTOFF) ) {
      xbar += fsdb->fss[i]->seq_len;
      ybar += fsdb->fss[i]->score;
      j++;
    }
  }
  xbar /= j ;
  ybar /= j ;

  double ssxy = 0, ssxx = 0 ;
  for ( i = 0; i < fsdb->num_fss; i++ ) {
    if ( fsdb->fss[i]->unique_best &&
	(fsdb->fss[i]->score >= FIRST_ROUND_SCORE_CUTOFF) ) {
      ssxy += (fsdb->fss[i]->seq_len - xbar) * (fsdb->fss[i]->score - ybar) ;
      ssxx += (fsdb->fss[i]->seq_len - xbar) * (fsdb->fss[i]->seq_len - xbar) ;
    }
  }
  *slope = ssxy / ssxx ;
  *intercept = ybar - *slope * xbar ;

  if (DEBUG) {
    FILE* LVSLOG = fileOpen( "LENvSCORE.dat", "w" );
    fprintf( LVSLOG,
	"# Just calculated length-score best-fit line:\n" );
    fprintf( LVSLOG,
	"# score = %0.4f + (length x %0.4f)\n",
	*intercept, *slope );
    for ( i = 0; i < fsdb->num_fss; i++ ) {
      if ( fsdb->fss[i]->unique_best &&
	  (fsdb->fss[i]->score >= FIRST_ROUND_SCORE_CUTOFF) ) {
	fprintf( LVSLOG, "%d\t%d\n", fsdb->fss[i]->seq_len, fsdb->fss[i]->score );
      }
    }
    fclose( LVSLOG );
  }
}

/* set_uniq_in_fsdb
   Args: (1) FSDB fsdb - has fss field SORTED!
   Returns: void
   Goes through each sequence and the sets the unique_best flag
   to true for the first of each kind (same as, ae, and rc) and
   sets unique_best to false for all others
*/
void set_uniq_in_fsdb( FSDB fsdb ) {
  int i, curr_rc, curr_as, curr_ae;
  FragSeqP fs;
  fs = fsdb->fss[0];
  /* initialize */
  curr_rc = fs->rc;
  curr_as = fs->as;
  curr_ae = fs->ae;
  fs->unique_best = 1;
  for ( i = 1; i < fsdb->num_fss; i++ ) {
    fs = fsdb->fss[i];

    /* If new guy is same as last guy, on strand, start, and end,
       he's redundant (not unique) */
    if ( (fs->rc == curr_rc) &&
	 (fs->as == curr_as) &&
	 (fs->ae == curr_ae) ) {
      fs->unique_best = 0;
    }

    else {
      /* forward strand */
      if ( fs->rc == 0 ) {
	/* Can still be redundant if it's untrimmed */
	if ( fs->as == curr_as ) {
	  if ( fs->trimmed ) {
	    /* strand and start match, end does not and it's trimmed:
	       therefore, it's a unique best */
	    fs->unique_best = 1;
	    curr_ae = fs->ae;
	  }
	  else {
	    fs->unique_best = 0;
	  }
	}
	else {
	  curr_rc = fs->rc;
	  curr_as = fs->as;
	  curr_ae = fs->ae;
	  fs->unique_best = 1;
	}
      }

      /* reverse strand */
      else {
	if ( fs->ae == curr_ae ) {
	  if ( fs->trimmed ) {
	    /* strand and end (beginning of rc molecule, dummy) match
	       but start (end, really) and it's trimmed, therefore
	       it's a unique best */
	    fs->unique_best = 1;
	    curr_as = fs->as;
	  }
	  else {
	    fs->unique_best = 0;
	  }
	}

	else {
	  curr_rc = fs->rc;
	  curr_as = fs->as;
	  curr_ae = fs->ae;
	  fs->unique_best = 1;
	}
      }
    }
  }
}


/* pop_smp_from_FSDB
   Args: (1) FSDB fsdb with valid data
         (2) Depth of PSSM matrices
   Returns: void
   Goes through all the AlnSeqs in the fsdb->fss array. Follows
   the front_asp (and back_asp, if necessary) pointer to populate
   the smp field of all AlnSeqs with the correct code for what
   depth in the PSSM matrix to use for consensus calling */
void pop_smp_from_FSDB( FSDB fsdb, int depth ) {
  int i, aln_seq_pos, front_seq_len, back_seq_len,
    distance_from_front,
    distance_from_back,
    aln_seq_len;
  int act_seq_pos = 0;
  AlnSeqP front_asp, back_asp;

  for ( i = 0; i < fsdb->num_fss; i++ ) {
    front_asp = fsdb->fss[i]->front_asp;
    back_asp  = fsdb->fss[i]->back_asp;
    act_seq_pos = 0;
    front_seq_len = asp_len( front_asp );
    if ( back_asp != NULL ) {
      back_seq_len  = asp_len( back_asp );
    }
    else {
      back_seq_len = 0;
    }

    /* First, fill in the front_asp->smp array */
    aln_seq_len = front_asp->end - front_asp->start + 1;
    for( aln_seq_pos = 0; aln_seq_pos < aln_seq_len; aln_seq_pos++ ) {
      if ( front_asp->ins[aln_seq_pos] != NULL ) {
	act_seq_pos += strlen( front_asp->ins[aln_seq_pos] );
      }
      distance_from_front = act_seq_pos;
      distance_from_back  = (front_seq_len + back_seq_len) -
	act_seq_pos - 1;

      if ( distance_from_front <= depth ) {
	front_asp->smp[aln_seq_pos] = ('A'+distance_from_front);
      }
      else {
	if ( distance_from_back < depth ) {
	  front_asp->smp[aln_seq_pos] = ('A'+(depth*2)-distance_from_back);
	}
	else {
	  front_asp->smp[aln_seq_pos] = ('A' + depth);
	}
      }

      if ( front_asp->seq[aln_seq_pos] != '-' ) {
	act_seq_pos++;
      }
    }
    front_asp->smp[aln_seq_pos] = '\0';

    /* Then, fill in the back_asp->smp array */
    if ( back_asp != NULL ) {
      aln_seq_len = back_asp->end - back_asp->start + 1;
      for( aln_seq_pos = 0; aln_seq_pos < aln_seq_len; aln_seq_pos++ ) {
	if ( back_asp->ins[aln_seq_pos] != NULL ) {
	  act_seq_pos += strlen( back_asp->ins[aln_seq_pos] );
	}
	distance_from_front = (front_seq_len + act_seq_pos);
	distance_from_back = (front_seq_len + back_seq_len) -
	  act_seq_pos - 1;
	if ( distance_from_front <= depth ) {
	  back_asp->smp[aln_seq_pos] = ('A' + distance_from_front);
	}
	else {
	  if ( distance_from_back < depth ) {
	    back_asp->smp[aln_seq_pos] = ('A'+(depth*2)-distance_from_back);
	  }
	  else {
	    back_asp->smp[aln_seq_pos] = ('A' + depth);
	  }
	}

	if ( back_asp->seq[aln_seq_pos] != '-' ) {
	  act_seq_pos++;
	}
      }
      back_asp->smp[aln_seq_pos] = '\0';
    }
  }
}


/* add_fs2fsdb
   Args: (1) FragSeqP fs - pointer to a fully valid FragSeq
         (2) FSdb fsdb - database to add this FragSeq to
   Returns: 1 if success; 0 if failure (not enough memories)
   Adds the FragSeq pointed to by fs to the fsdb database,
   growing it if necessary */
int add_fs2fsdb( FragSeqP fs, FSDB fsdb ) {
  FragSeqP next_fs;

  /* First, check if fsdb need to grow */
  if ( fsdb->num_fss == fsdb->size ) {
    if ( grow_FSDB( fsdb ) == 0 ) {
      return 0;
    }
  }

  /* Get a pointer to the next available FragSeq */
  next_fs = fsdb->fss[fsdb->num_fss];

  /* Copy over the input fs into next_fs */
  strcpy( next_fs->id, fs->id );
  strcpy( next_fs->desc, fs->desc );
  strcpy( next_fs->seq, fs->seq );
  next_fs->trim_point = fs->trim_point;
  next_fs->trimmed    = fs->trimmed;
  next_fs->seq_len    = fs->seq_len;
  next_fs->rc         = fs->rc;
  next_fs->as         = fs->as;
  next_fs->ae         = fs->ae;
  next_fs->score      = fs->score;
  next_fs->front_asp  = fs->front_asp;
  next_fs->back_asp   = fs->back_asp;
  next_fs->unique_best = fs->unique_best;

  /* Bump up the num_fss */
  fsdb->num_fss += 1;
  return 1;
}

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
int grow_FSDB( FSDB fsdb ) {
  int i, j, new_size;
  FragSeqP first_seq;
  FragSeqP* new_fss;

  new_size = fsdb->size * 2;

  /* DEBUG INFO */
  if ( DEBUG ) {
    fprintf( stderr, "Growing fsdb from %d to %d\n",
	     fsdb->size, new_size );
  }

  /* Allocate another chunck of memories as big as the
     one it has now, doubling its size */
  first_seq = (FragSeqP)save_malloc(fsdb->size *
			       sizeof(FragSeq));
  if ( first_seq == NULL ) {
    return 0;
  }

  /* Now, allocate the *new* array of pointers for fsdb->fss
     But, assign this to new_fss for now because we need to
     keep fsdb->fss so we can copy over the pointers it already
     has!
  */
  new_fss = (FragSeqP*)save_malloc(new_size * sizeof(FragSeqP));
  if ( new_fss == NULL ) {
    return 0;
  }

  /* Point the pointers to the pointees */
  for( i = 0; i < fsdb->size; i++ ) {
    new_fss[i] = fsdb->fss[i];
  }
  j = 0;
  for( i = fsdb->size; i < new_size; i++ ) {
    new_fss[i] = &first_seq[j++];
  }

  /* Now, free the old fsdb->fss and slot in the new one */
  free( fsdb->fss );
  fsdb->fss  = new_fss;
  fsdb->size = new_size;
  return 1;
}

/* init_FSDB
   Arguments: void
   Returns: FSDB (pointer to struct fragseqdb) / NULL if
    not enough memories
   Used for initializing a new database of FragSeqs. Allocates
   enough memoery for INIT_NUM_ALN_SEQS of these
*/
FSDB init_FSDB ( void ) {
  int i;
  FSDB fsdb;
  FragSeqP first_seq;

  /* First, allocate the memories */
  fsdb = (FSDB)save_malloc(sizeof(FragSeqDB));
  if ( fsdb == NULL ) {
    return NULL;
  }

  first_seq = (FragSeqP)save_malloc(INIT_NUM_ALN_SEQS *
			       sizeof(FragSeq));
  if ( first_seq == NULL ) {
    return NULL;
  }

  fsdb->fss = (FragSeqP*)save_malloc(INIT_NUM_ALN_SEQS *
				sizeof( FragSeqP ));
  if ( fsdb->fss == NULL ) {
    return NULL;
  }

  for ( i = 0; i < INIT_NUM_ALN_SEQS; i++ ) {
    fsdb->fss[i] = &first_seq[i];
  }

  fsdb->size = INIT_NUM_ALN_SEQS;
  fsdb->num_fss = 0;

  return fsdb;
}

