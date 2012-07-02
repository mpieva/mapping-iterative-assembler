/* $Id$ */
#include "mia.h"

/*
    void * save_malloc(size_t size){
        void * tmp = malloc(size);
        if (tmp == NULL)
            fprintf( stderr, "Some memory allokation failed. Exiting");

        return tmp;
    }
*/

/* c2rcc
   Args:(1) int c - the coordinate (0-indexed)
        (2) int len - the length of the sequence on which c is
   Returns: an int - the coordinate on the reverse complement
   This function returns the corresponding coordinate for any
   position on a sequence on the reverse complement of that
   sequence. It also works if the sequence and its reverse
   complement have been "wrapped", i.e., a bit of sequence
   from the beginning has been added to the end so long as
   the input len is of the *ORIGINAL* sequence the was then
   wrapped
*/
int c2rcc ( const int c, const int len ) {
  int mc;
  mc = c%len; /* find c mod len */
  return len - mc - 1;
}

/* init_sized_map_alignment
   Args: MapAlignmentP maln - source maln to get size from
   Returns: pointer to a fresh MapAligment for holding the
   culled results from the source maln. This culled guy
   gets the same ref and an AlnSeqArray big enough for the
   results, but no actual memory is malloced for new AlnSeq's
   Instead, we'll just point the source maln guys to this 
   one if they are unique, which is determined elsewhere
*/
MapAlignmentP init_culled_map_alignment( MapAlignmentP src_maln ) {
  MapAlignmentP culled_maln;
  culled_maln = (MapAlignmentP)save_malloc(sizeof(MapAlignment));
  if ( culled_maln == NULL ) {
    return NULL;
  }
  culled_maln->ref = src_maln->ref;
  culled_maln->AlnSeqArray = (AlnSeqP*)save_malloc(src_maln->num_aln_seqs *
					      sizeof(AlnSeqP));
  if ( culled_maln->AlnSeqArray == NULL ) {
    return NULL;
  }
  culled_maln->num_aln_seqs = 0;
  culled_maln->size = src_maln->num_aln_seqs;
  culled_maln->cons_code = src_maln->cons_code;
  culled_maln->distant_ref = src_maln->distant_ref;
  return culled_maln;
}

/* find_alignable_len
   Args: (1) FragSeqP fs - with value info in as, ae and seq_len fields
         (2) RefSeqP ref - with valid info in the sequence
   Returns: int with the alignable sequence length of this sequence
   in this FragSeqP. That is defined as the length of this sequence minus
   any part that overlaps positions that are "N" in the RefSeq. This
   number is not allowed to be less that MIN_ALIGNABLE_LEN to avoid
   having sequence with very little or no alignable sequence.
*/
int find_alignable_len( FragSeqP fs, RefSeqP ref ) {
  int alignable_len;
  size_t i;
  alignable_len = fs->seq_len;

  for ( i = fs->as; i < fs->ae; i++ ) {
    if ( ref->seq[i] == 'N' ) {
      alignable_len--;
    }
  }

  if ( alignable_len < MIN_ALIGNABLE_LEN ) {
    alignable_len = MIN_ALIGNABLE_LEN;
  }

  return alignable_len;
}

inline int valid_base( char b ) {
  if ( (b == 'A') ||
       (b == 'C') ||
       (b == 'G') ||
       (b == 'T') ) {
    return 1;
  }
  else {
    return 0;
  }
}

/* NEED TO CHANGE THE WAY Q-SCORES ARE ADDED AND SUBTRACTED! -
   Done (Ed 5 Oct 2009)
*/
/* collapse_fs
   Args: (1) FragSeqP dest_fs - destination FragSeqP to be modified
         (2) FragSeqP fs_to_add - FragSeqP to be melded into dest_fs
	                          and then marked for removal
   Returns: void
   This function is called when two FragSeqP's are identified that
   are redundant. They are collapsed together by going through, base
   by base, and checking for the highest quality score. If the bases
   match, quality scores are added. If they differ, the base with 
   the higher q-score is set in the dest_fs and the quality scores
   are subtracted. Non-base characters are ignored since they give
   no information.
   If the sequences are of different length, collapsing is not done.
   Quality scores are capped at 255
 */
void collapse_fs( FragSeqP dest_fs, FragSeqP fs_to_add ) {
  size_t i, seq_len;
  int test_qual;
  char dest_base, add_base;
  /* Don't collapse if sequences are not the same length! */
  if ( dest_fs->seq_len != fs_to_add->seq_len ) {
    return;
  }
  seq_len = dest_fs->seq_len;
  for( i = 0; i < seq_len; i++ ) {
    dest_base = dest_fs->seq[i];
    add_base  = fs_to_add->seq[i];
    test_qual = 0;
    /* If either base is not valid, just use the valid one
       verbatem */
    if ( !valid_base( add_base ) ) {
      1;
    }
    else {
      if ( !valid_base( dest_base ) ) {
	dest_fs->seq[i]  = add_base;
	dest_fs->qual[i] = dest_fs->qual[i];
      }

      /* Both are valid bases, let's collapse! */
      else {
	if (dest_base == add_base) {
	  test_qual = dest_fs->qual[i] + (fs_to_add->qual[i] - QUAL_ASCII_OFFSET);
	  //dest_fs->qual[i] += (fs_to_add->qual[i] - QUAL_ASCII_OFFSET);
	}
	else {      /* Disagreement */
	  
	  /* Better quality from the destination? */
	  if ( dest_fs->qual[i] >= fs_to_add->qual[i] ) {
	    test_qual = dest_fs->qual[i] - 
	      (fs_to_add->qual[i] - QUAL_ASCII_OFFSET);
	    //dest_fs->qual[i] -= (fs_to_add->qual[i] - QUAL_ASCII_OFFSET);
	  }
	  /* Better quality from the new guy, switch-a-roo */
	  else {
	    dest_fs->seq[i] = add_base;
	    test_qual = fs_to_add->qual[i] - 
	      (dest_fs->qual[i] - QUAL_ASCII_OFFSET);
	    //	    dest_fs->qual[i] = fs_to_add->qual[i] - 
	    //	      (dest_fs->qual[i] - QUAL_ASCII_OFFSET);
	  }
	}
      }
      /* Range check */
      if ( test_qual >= 127 ) {
	dest_fs->qual[i] = '~';
      }
      if ( test_qual <= 33 ) {
	dest_fs->qual[i] = '!';
      }
      /*      if ( dest_fs->qual[i] >= 127 ) {
	dest_fs->qual[i] = 127;
	}*/
    }
  }
  dest_fs->num_inputs++;
  fs_to_add->num_inputs = 0;
}

/* collapse_FSDB
   Args: (1) FSDB fsdb
   Returns: void
   Goes through each FragSeq pointed to by fsdb->fss. They must be sorted
   and have their unique_best field populated. For all that pass
   the score cutoff filter and match in the rc, as, and ae fields, 
   collapses them into a single sequence. Updates their quality scores
   to be aggregates.
*/
void collapse_FSDB( FSDB fsdb, int Hard_cut, 
		    int SCORE_CUT_SET, double s, double n ) {
  FragSeqP fs, current_collapsing_fs;
  size_t i, j;
  double slope_def = 100.0;
  double slope, intercept;
  double min_score_for_len;
   /* Using user-defined slope and intercept */
  if ( SCORE_CUT_SET ) {
    slope = s;
    intercept = n;
  }
  else {
    /* Find best-fit for length-score relationship */
    find_fsdb_score_cut( fsdb, &slope, &intercept );
  }

  /* Make sure these are sensible values */
  if ( slope <= 0 ) {
    slope = slope_def;
  }
 
  for ( i = 0; i < fsdb->num_fss; i++ ) {
    fs = fsdb->fss[i];
    if ( SCORE_CUT_SET ) {
      min_score_for_len = (double)(intercept + (slope * fs->seq_len));
    }
    else {
      min_score_for_len = (double)(intercept + (slope * fs->seq_len));
    }
    if ( Hard_cut > 0 ) {
      min_score_for_len = Hard_cut;
    }

    if ( fs->unique_best &&
	 (fs->score >= min_score_for_len) ) {
      current_collapsing_fs = fs;
      i++;
      while( (i < fsdb->num_fss) &&
	     (fsdb->fss[i]->unique_best == 0) ) {
	fs = fsdb->fss[i];
	if ( fs->score >= min_score_for_len ) {
	  /* Collapsing puts these two together, increments
	     the num_inputs field in the current_collapsing_fs,
	     and sets the num_inputs field to 0 in the fs to
	     mark it for removal */
	  collapse_fs( current_collapsing_fs, fs );
	}
	i++;
      }
      i--;
    }
  }

  /* Time to take out the trash */
  j = 0; // j will be the next available slot
  for ( i = 0; i < fsdb->num_fss; i++ ) {
    if ( fsdb->fss[i]->num_inputs > 0 ) {
      fsdb->fss[j] = fsdb->fss[i];
      j++;
    }
  }
  fsdb->num_fss = j;
}

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
			  int SCORE_CUT_SET, double s, double n ) {
  int i, j, ref_gaps, new_ref_gaps, culled_nas, alignable_len;
  FragSeqP fs;
  AlnSeqP aln_seq;
  char* ins_seq;
  double slope_def = 100.0;
  double slope, intercept; /* Parameters for score-culling */
  double min_score_for_len;

  /* Using user-defined slope and intercept */
  if ( SCORE_CUT_SET ) {
    slope = s;
    intercept = n;
  }
  else {
    /* Find best-fit for length-score relationship */
    find_fsdb_score_cut( fsdb, &slope, &intercept );
  }

  /* Make sure these are sensible values */
  if ( slope <= 0 ) {
    slope = slope_def;
  }

  /* Print these, if DEBUG */
  if ( DEBUG ) {
    fprintf( stderr, "Culling maln using slope=%0.4f, intercept=%0.4f, percent of avg.=%d\n",
	     slope, intercept, SCORE_CUTOFF_BUFFER );
  }

  culled_nas = 0;
  for ( i = 0; i < fsdb->num_fss; i++ ) {
    fs = fsdb->fss[i];
    if ( SCORE_CUT_SET ) {
      min_score_for_len = (double)(intercept + (slope * fs->seq_len));
    }
    else {
      /* If we're using a distant reference genome, we may have sequence
	 that overlaps many N positions. These should not be counted
	 againt the length of this sequence */
      if ( culled_maln->distant_ref ) {
	alignable_len = find_alignable_len( fs, culled_maln->ref );
	min_score_for_len =  (double)(intercept + (slope * alignable_len));
	  
      }
      else {
	min_score_for_len =  (double)(intercept + (slope * fs->seq_len));
      }
    }
    if ( Hard_cut > 0 ) {
      min_score_for_len = Hard_cut;
    }
    if ( fs->unique_best && 
	 (fs->score >= min_score_for_len) ) {
      culled_maln->AlnSeqArray[culled_nas++] = fs->front_asp;
      if ( fs->back_asp != NULL ) {
	culled_maln->AlnSeqArray[culled_nas++] = fs->back_asp;
      }
    }
  }
  culled_maln->num_aln_seqs = culled_nas;

  /* Now, reset the culled_maln->ref->gaps array since some of these
     may have come from sequences that are now removed */
  for( i = 0; i < culled_maln->ref->seq_len; i++ ) {
    ref_gaps = culled_maln->ref->gaps[i];
    if ( ref_gaps > 0 ) {
      new_ref_gaps = 0;
      for( j = 0; j < culled_maln->num_aln_seqs; j++ ) {
	aln_seq = culled_maln->AlnSeqArray[j];
	if ( (aln_seq->start <  i) &&
	     (aln_seq->end   >= i) ) {
	  /* Does it have some actual inserted sequence? */
	  ins_seq = aln_seq->ins[i - aln_seq->start];
	  if ( (ins_seq != NULL) &&
	       (strlen( ins_seq ) > new_ref_gaps ) ) {
	    new_ref_gaps = strlen( ins_seq );
	  }
	}
      }
      culled_maln->ref->gaps[i] = new_ref_gaps;
    }
  }
  return;
}

/* asp_len
   Args: (1) AlnSeqP asp - pointer to an AlnSeq
   Returns (1) int - total length of sequence 
   This function finds the total length of the sequence in this
   aligned sequence fragment. This is the sum of the sequence
   in the asp->seq field and all of the inserted sequence (if any)
   in the asp->ins array
*/
int asp_len( AlnSeqP asp ) {
  int i;
  int aln_seq_len = 0;
  int tot_seq_len = 0;
  aln_seq_len = (asp->end - asp->start + 1);
  tot_seq_len = (asp->end - asp->start + 1);
  for( i = 0; i < aln_seq_len; i++ ) {
    if ( asp->ins[i] != NULL ) {
      tot_seq_len += strlen( asp->ins[i] );
    }
  }
  return tot_seq_len;
}


/* consensus_assembly_string
   Takes a maln object
   Generates the consensus sequence string using the aligned data
   within the maln according to the maln->cons_code and puts it
   in char* cons
   Returns char* pointer to consensus string
*/
char* consensus_assembly_string ( MapAlignmentP maln ) {

  int cons_pos, ref_pos, ref_gaps, j, num_gaps;
  char ins_cons[MAX_INS_LEN + 1];
  int  ins_cov[MAX_INS_LEN + 1];
  char cons_base;
  char* cons;
  AlnSeqP aln_seq; // pointer to currently considered aligned sequence
  BaseCountsP bcs;
  PSSMP psm;
  
  num_gaps = 0;

  for( j = 0; j < maln->ref->seq_len; j++ ) {
    num_gaps += maln->ref->gaps[j];
  }

  cons = (char*)save_malloc((maln->ref->seq_len + num_gaps +1) *
		       sizeof(char));
  if ( cons == NULL){
      fprintf( stderr, "Not enough memory for cons\n" );
      exit(1);
  }

  bcs = (BaseCountsP)save_malloc(sizeof(BaseCounts));
  if ( bcs == NULL ) {
    fprintf( stderr, "Not enough memories for bcs\n" );
    exit( 1 );
  }
  reset_base_counts(bcs);

  cons_pos = 0;
  ref_pos = 0;

  /* Go through each position of the reference sequence, as it
     currently is */
  for ( ref_pos = 0; ref_pos < maln->ref->seq_len; ref_pos++ ) {
    /* How many gaps preceed this position? */
    ref_gaps = maln->ref->gaps[ref_pos];
    if ( (ref_gaps > 0) &&
	 (ref_pos  > 0) ) {
      find_ins_cons( maln, ref_pos, ins_cons, ins_cov, 0 );
      for ( j = 0; j < ref_gaps; j++ ) {
	/* Consensus is a gap, i.e., nothing. So, do not write this
	   to the consensus assembly */
	if ( (ins_cons[j] == '-') ||
	     (ins_cons[j] == ' ') ) {
	  ;
	}
	else {
	  cons[cons_pos] = ins_cons[j];
	  cons_pos++;
	}
      }
    }

    /* Re-zero all base counts */
    reset_base_counts(bcs);

    /* Find all the aligned fragments that include this
       position and make a consensus from it */
    for( j = 0; j < maln->num_aln_seqs; j++ ) {
      aln_seq = maln->AlnSeqArray[j];

      /* Does this aligned fragment cover this position */
      if ( (aln_seq->start <= ref_pos) && // checked
	   (aln_seq->end   >= ref_pos) ) {

	if ( aln_seq->revcom ) {
	  psm = maln->rpsm;
	}
	else {
	  psm = maln->fpsm;
	}

	add_base( aln_seq->seq[ref_pos - aln_seq->start], 
		  bcs, psm,
		  aln_seq->smp[ref_pos - aln_seq->start] );
      }
    }
    cons_base = find_consensus( bcs, maln->cons_code );
    if ( (cons_base == '-') ||
	 (cons_base == ' ') ) {
      ;
    }
    else {
      cons[cons_pos] = cons_base; 
      cons_pos++;
    }
  }
  cons[cons_pos] = '\0';
  free( bcs );
  return cons;
}

/* Takes a pointer to an Alignment that has valid values
   in its dynamic programming matrix and valid values
   for a->aer and a->aec (ending row and column).
   Tracks back to the beginning of the alignment and
   adds valid values from a->abr and a->abc
   Returns nothing
*/
void find_align_begin( AlignmentP a ) {
  int row, col;
  row = a->aer;
  col = a->aec;
  
  while( (a->m->mat[row][col].trace != col) &&
	 (a->m->mat[row][col].trace != -row) ) {
    if ( a->m->mat[row][col].trace == 0 ) {
      row--;
      col--;
    }
    else {
      if ( a->m->mat[row][col].trace < 0 ) {
	row = -(a->m->mat[row][col].trace);
	col--;
      }
      else {
	col = a->m->mat[row][col].trace;
	row--;
      }
    }
  }

  a->abc = col;
  a->abr = row;
}




void make_ref_upper( RefSeqP ref ) {
  size_t pos;
  for( pos = 0; pos < ref->wrap_seq_len; pos++ ) {
    ref->seq[pos]   = toupper( ref->seq[pos] );
    ref->rcseq[pos] = toupper( ref->rcseq[pos] );
  }
}

/* Takes a pointer to a RefSeq that has a valid 
   sequence in it. Adds INIT_ALN_SEQ_LEN sequence
   from the beginning to the end so that any
   sequence fragment aligned to it will have a
   valid chance to align, despite the circularity
   of the sequence. Returns TRUE if this operation
   went smoothly, FALSE otherwise */
int add_ref_wrap( RefSeqP ref ) {
  int wrap_len; // amount of sequence to wrap around
  
  /* Can't wrap more than is there! */
  if ( ref->seq_len < INIT_ALN_SEQ_LEN ) {
    wrap_len = ref->seq_len;
  }
  else {
    wrap_len = INIT_ALN_SEQ_LEN;
  }

  /* Grow the ref->seq and ref->rcseq if they're not already big enough */
  while ( (ref->seq_len + wrap_len) >= ref->size ) {
    ref->seq = grow_seq( ref->seq, ref->size );
    if ( ref->rcseq != NULL ) {
      ref->rcseq = grow_seq( ref->rcseq, ref->size );
    }
    ref->size = ref->size * 2;
  }

  /* Now copy wrap_len characters from the beginning of
     ref to the end of ref */
  strncat( &ref->seq[ref->seq_len], ref->seq, wrap_len );
  if ( ref->rcseq != NULL ) {
    strncat( &ref->rcseq[ref->seq_len], ref->rcseq, wrap_len );
  }

  /* set the wrap_seq_len */
  ref->wrap_seq_len = wrap_len + ref->seq_len;

  /* set the circular flag */
  ref->circular = 1;
}

/* init_dpm
   Args: (1) size1 - the number of rows (fragment sequence)
         (2) size2 - the number of columns (referense sequence)
   Returns: DPMP => pointer to a dynamic programming matrix
   with memory properly allocated
*/
DPMP init_dpm( int size1, int size2 ) {
  DPMP m;
  DPEP elements;
  int i;

  m = (DPMP)save_malloc(sizeof(Mat));
  m->rows = size1;
  m->cols = size2;

  /* Allocate the elements */
  elements = (DPEP)save_malloc(m->rows * m->cols * sizeof(DPE));
  if ( elements == NULL ) {
    return NULL;
  }

  /* Allocate the rows */
  m->mat = (DPEP*)save_malloc(m->rows * sizeof(DPEP));
  if ( m->mat == NULL ) {
    free(elements);
    return NULL;
  }

  /* Assign the rows the proper values */
  for ( i = 0; i < m->rows; i++ ) {
    m->mat[i] = &elements[m->cols * i];
  }
  return m;
}

void free_dpm( DPMP m ) {
  if( m ) {
    free( m->mat[0] ) ;
    free( m->mat ) ;
  }
  free( m ) ;
}


/* Takes a pointer to an Alignment
   that has valid sequence, length, submat, and sg data
   Does dynamic programming, filling in values in the 
   a->m dynamic programming matrix
   Returns nothing */
void dyn_prog( AlignmentP a ) {
  int row, 
    col, 
    gap_col_score, 
    gap_row_score,
    diag_score,
    start_new_score,
    hp_disc_gap_col_score,
    hp_disc_gap_row_score,
    sm_depth;
  size_t i;
  int HIM = (INT_MIN / 2); // half of the minimum int; this
  //is useful to avoid underflow from subtracting from the
  // smallest possible int
  int row_sm[5]; // row substitution matrix
  
  /* Initialize */
  row = 0;
  col = 0;
  hp_disc_gap_col_score = HIM;
  hp_disc_gap_row_score = HIM;

  /* Set up the substitution matrix scores for this base for
     this first row */
  for( i = 0; i <= 4; i++ ) {
    row_sm[i] = a->submat->sm[0][i][a->s2c[0]];
  }

  // First row, no penalty whether sg or not
  for( col = 0; col < a->len1; col++ ) {
    if ( a->align_mask[col] ) {
      a->m->mat[row][col].score =
	row_sm[a->s1c[col]];

      //      a->m->mat[row][col].score = 
      //sub_mat_score(a->s1c[col], a->s2c[row],
      //	      a->submat->sm, row, a->len2);
    }
    else {
      a->m->mat[row][col].score = HIM;
    }
    a->m->mat[row][col].trace = 0; // any alignment
    //traced back this far must start here, just stick
    //in a 0 for the heck of it
    a->best_gap_row[col] = 0;
  }

  // Subsequent rows, 
  for ( row = 1; row < a->len2; row++ ) {
    
    /* Set up the substitution matrix scores for the base 
       for this row */
    sm_depth = find_sm_depth( row, a->len2 );
    for( i = 0; i <= 4; i++ ) {
      row_sm[i] = a->submat->sm[sm_depth][i][a->s2c[row]];
    }
    /* First column, special case - no gapping
       in either direction and pay a penalty for
       unaligned bases if sg */
    col = 0;
    /*
      assert( a && a->m && a->m->mat && a->s1c && a->s2c && a->submat ) ;
      assert( 0 <= row && row < a->len2 ) ;
      assert( 0 <= col && col < a->len1 ) ;
    */
    if ( a->align_mask[col] ) {
      /* a->m->mat[row][col].score = 
	sub_mat_score(a->s1c[col], a->s2c[row],
	a->submat->sm, row, a->len2); */
      a->m->mat[row][col].score =
	row_sm[a->s1c[col]];
      // pay penalty at col 0 if sg
      if ( a->sg5 ) {
	a->m->mat[row][col].score 
	  -= (GOP + (GEP * (row+1)));
      }
    }

    else {
      a->m->mat[row][col].score = HIM;
    }

    a->m->mat[row][col].trace = 0;

    // Subsequent columns
    a->best_gap_col = 0;
    for( col = 1; col < a->len1; col++ ) {
      if ( a->align_mask[col] ) {
	/*	a->m->mat[row][col].score = 
	  sub_mat_score(a->s1c[col], a->s2c[row],
	  a->submat->sm, row, a->len2); */
	a->m->mat[row][col].score =
	  row_sm[a->s1c[col]];
	
	/* update best_gap_col by comparing new gap
	   option to previous best, if we're far 
	   enough in to gap columns
	*/
	if ( col >= 2 ) {
	  if ( (a->m->mat[row-1][col-2].score - (GOP + GEP)) >
	       (a->m->mat[row-1][a->best_gap_col].score - 
		(GOP + (GEP * (col-(a->best_gap_col)-1)))) ) {
	    a->best_gap_col = col - 2;
	  }
	  gap_col_score = 
	    ( a->m->mat[row-1][a->best_gap_col].score -
	      (GOP + (GEP * (col - a->best_gap_col - 1))) );
	}
	else {
	  gap_col_score = HIM;
	}

	/* update best_gap_row by comparing new gap
	   option to previous best, if we're far enough
	   down to gap rows
	*/
	if ( row >= 2 ) {
	  if ( (a->m->mat[row-2][col-1].score - (GOP + GEP)) > 
	       (a->m->mat[a->best_gap_row[col-1]][col-1].score -
		(GOP + (GEP * (row-(a->best_gap_row[col-1])-1)))) ) {
	    a->best_gap_row[col-1] = row - 2;
	  }
	  gap_row_score = 
	    ( a->m->mat[a->best_gap_row[col-1]][col-1].score -
	      (GOP + (GEP * (row-(a->best_gap_row[col-1])-1))) );
	}
	else {
	  gap_row_score = HIM;
	}
	
	/* Find diagonal score */
	diag_score = a->m->mat[row-1][col-1].score;
	
	/* Check if starting a new alignment here is the best
	   option. If it's a->sg5 TRUE, then we have to pay
	   the penalty for unaligned beginning sequence. If
	   not, no penalty and we simply start from zero */
	start_new_score = 0;
	if ( a->sg5 ) {
	  start_new_score -= (GOP + (GEP * (row+1)));
	}
	
	/* Calculate special homopolymer discount? */
	if ( a->hp ) {
	  hp_disc_gap_col_score = HIM;
	  hp_disc_gap_row_score = HIM;
	  if ( a->seq1[col] == a->seq2[row] ) { // must be the same base
	    if ( (a->hprs[row] == row) && // seq1 hp starts here
		 (a->hpcs[col] != col) && // seq2 hp starts before here
		 (a->hpcs[col] > 0) // can't gap outside of seq1!
		 ) {
	      hp_disc_gap_col_score = 
		(a->m->mat[row-1][(a->hpcs[col]-1)].score -
		 hp_discount_penalty( (col - a->hpcs[col]),
				      a->hpcl[col], a->hprl[row] ));
	    }
	    if ( (a->hpcs[col] == col) && // seq2 hp starts here
		 (a->hprs[row] != row) && // seq1 hp starts before here
		 (a->hprs[row] > 0) ) { // can't gap outside of seq2!
	      hp_disc_gap_row_score = 
		(a->m->mat[(a->hprs[row]-1)][col-1].score -
		 hp_discount_penalty( (col - a->hpcs[col]),
				      a->hpcl[col], a->hprl[row] ));
	    }
	  }
	}
	
	/* Now best options are on the table, pick the
	   best of the best */
	/* Best option is starting a new alignment */
	if ( (start_new_score > diag_score) &&
	     (start_new_score > gap_col_score) &&
	     (start_new_score > gap_row_score) &&
	     (start_new_score > hp_disc_gap_col_score) &&
	     (start_new_score > hp_disc_gap_row_score)
	     ) {
	  a->m->mat[row][col].trace = col; /* Mark this as the beginning */
	  a->m->mat[row][col].score = start_new_score;
	}
	
	else {
	  /* Best option is continuing on the diagonal */
	  if ( (diag_score >= gap_col_score) &&
	       (diag_score >= gap_row_score) &&
	       (diag_score >= hp_disc_gap_col_score) &&
	       (diag_score >= hp_disc_gap_row_score)
	       ) {
	    a->m->mat[row][col].trace = 0;
	    a->m->mat[row][col].score += diag_score;
	  }
	  
	  else {
	    /* Is best option gapping back columns? */
	    if ( (gap_col_score >= gap_row_score) &&
		 (gap_col_score >= hp_disc_gap_col_score) &&
		 (gap_col_score >= hp_disc_gap_row_score)
		 ) {
	      a->m->mat[row][col].score += gap_col_score;
	      a->m->mat[row][col].trace = a->best_gap_col;
	    }
	    
	    else {
	      if ( (gap_row_score >= hp_disc_gap_col_score) &&
		   (gap_row_score >= hp_disc_gap_row_score) ) {
		/* Best option must be gapping up rows */
		a->m->mat[row][col].score += gap_row_score;
		a->m->mat[row][col].trace = 
		  -(a->best_gap_row[col-1]);
	      }
	      else {
		if ( hp_disc_gap_col_score >= hp_disc_gap_row_score ) {
		  /* Best option is homopolymer discounted gapping
		     of columns */
		  a->m->mat[row][col].score += hp_disc_gap_col_score;
		  a->m->mat[row][col].trace = a->hpcs[col] - 1;
		}
		else {
		  /* Best option is homopolymer discounted gapping
		     of rows */
		  a->m->mat[row][col].score += hp_disc_gap_row_score;
		  a->m->mat[row][col].trace = -(a->hprs[row] - 1);
		}
	      }
	    }
	  }
	}
      }
      else {
	a->m->mat[row][col].score = HIM;
	a->m->mat[row][col].trace = 0;
      }
    }

    /* Now, end of row, so pay penalty for unaligned
       seq1, if sg3 alignment */
    if ( (a->sg3) &&
	 (a->len1 > (row+1)) ) {
      a->m->mat[row][col].score -= 
	GOP + ((a->len1 - row - 1) * GEP);
    }
  }
}

/* size1 is length of fragment
   size2 is length of reference (wrapped if necessary) + INIT_ALN_SEQ_LEN
   rc is boolean to seay if its reverse complement
   hp_special is boolean to say if homopolymer special gap costs are to be used
*/
AlignmentP init_alignment( int size1, int size2, 
			   int rc, int hp_special ) {
  AlignmentP al = (AlignmentP)save_malloc(sizeof(Alignment));
  al->gop = GOP;
  al->gep = GEP;
  al->hp  = hp_special;
  al->m   = init_dpm( size1, size2 );
  if ( al->m == NULL ) {
    return NULL;
  }

  /* Allocate memories for the alignment masks */
  al->align_mask = 
    (unsigned char*)save_malloc(size2 * 
			   sizeof(unsigned char));
  if ( al->align_mask == NULL ) {
    return NULL;
  }
  /* Set it up to be all unmasked by default */
  memset( al->align_mask, 1, size2 );

  al->s1c = (short int*)save_malloc(size2 * sizeof(short int));
  al->best_gap_row = (int*)save_malloc(size2 * sizeof(int));
  al->sg5 = 0; // initialize to local alignment
  al->sg3 = 0; // initialize to local alignment
  al->rc = rc; // set reverse complement boolean

  /* If user wants special hp gap discount, allocate
     memories for hpc1l, hpcs, hprl, and hprs */
  if ( al->hp ) {
    al->hpcl = (int*)save_malloc(size2*sizeof(int));
    al->hpcs = (int*)save_malloc(size2*sizeof(int));
    al->hprl = (int*)save_malloc(size1*sizeof(int));
    al->hprs = (int*)save_malloc(size1*sizeof(int));
  }
  else {
    al->hpcl = NULL;
    al->hpcs = NULL;
    al->hprl = NULL;
    al->hprs = NULL;
  }
  return al;
}

void free_alignment( AlignmentP al ) {
  if( al ) {
    free( al->hprs ) ;
    free( al->hprl ) ;
    free( al->hpcs ) ;
    free( al->hpcl ) ;
    free( al->best_gap_row ) ;
    free( al->s1c ) ;
    free_dpm( al->m ) ;
  }
  free( al ) ;
}


/* pop_s1c_in_a
   Args: (1) AlignmentP a - has a->seq1 and a->len1 set to 
   valid values
   Returns: void
   Populates the a->s1c array with code for quick lookup
   in submat
*/
void pop_s1c_in_a ( AlignmentP a ) {
  size_t i;
  int r_len;
  char b;
  r_len = a->len1;
  /*  if ( r_len < 1 ) {
    return;
    }*/

  for ( i = 0; i < r_len; i++ ) {
    b = a->seq1[i];
    switch( b ) {
    case 'A' :
      a->s1c[i] = 0;
      break;
    case 'C' :
      a->s1c[i] = 1;
      break;
    case 'G' :
      a->s1c[i] = 2;
      break;
    case 'T' :
      a->s1c[i] = 3;
      break;
    default:
      a->s1c[i] = 4;
    }
  }
}

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
int hp_discount_penalty ( int gap_len, int hplen1, int hplen2 ) {
  int penalty;
  penalty = ( GEP * gap_len );

  /* Trying function: 1/x in this table */
  switch ( hplen2 ) {
  case 1 :
    penalty += GOP;
    break;
  case 2 :
    penalty += GOP * 0.5;
    break;
  case 3 :
    penalty += GOP * 0.33;
    break;
  case 4 :
    penalty += GOP * 0.25;
    break;
  case 5 :
    penalty += GOP * 0.2;
    break;
  case 6 :
    penalty += GOP * 0.17;
    break;
  case 7 :
    penalty += GOP * 0.14;
    break;
  case 8 :
    penalty += GOP * 0.13;
    break;
  case 9 :
    penalty += GOP * 0.11;
    break;
  case 10 :
    penalty += GOP * 0.10;
    break;
  default :
    penalty += GOP * 0.10;
  }

  /* Trying function: 1/sqrt(x) in this table */
  /* In test of Neandertal versus human, didn't give high enough
     discount
  switch ( hplen2 ) {
  case 1 :
    penalty += GOP;
    break;
  case 2 :
    penalty += GOP * 0.71;
    break;
  case 3 :
    penalty += GOP * 0.57;
    break;
  case 4 :
    penalty += GOP * 0.50;
    break;
  case 5 :
    penalty += GOP * 0.45;
    break;
  case 6 :
    penalty += GOP * 0.41;
    break;
  case 7 :
    penalty += GOP * 0.38;
    break;
  case 8 :
    penalty += GOP * 0.35;
    break;
  case 9 :
    penalty += GOP * 0.33;
    break;
  case 10 :
    penalty += GOP * 0.32;
    break;
  default :
    penalty += GOP * 0.30;
  }
  */
  /* First try - penalty falls too fast */
  //  penalty = (GOP/hplen2) + (GEP * gap_len);

  return penalty;
}

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
void pop_hpl_and_hps ( const char* seq, int len, int* hpl, int* hps ) {
  char cur_base, last_base;
  int cur_pos, cur_hp_start, cur_hp_len, back_pos;
  if ( len < 1 ) {
    return;
  }

  /* Initialize */
  last_base    = seq[0];
  cur_hp_start = 0;
  cur_hp_len   = 1;
  hps[0] = 0;

  /* Go through each subsequenct position */
  for( cur_pos = 1; cur_pos < len; cur_pos++ ) {
    cur_base = seq[cur_pos];
    /* More of the same homopolymer? */
    if ( cur_base != last_base ) {
      /* New homopolymer, so back-fill the hpl */
      for( back_pos = cur_pos - 1; 
	   back_pos >= cur_hp_start; 
	   back_pos-- ) {
	hpl[back_pos] = cur_hp_len;
      }
      cur_hp_start = cur_pos;
      cur_hp_len   = 1;
    }
    else {
      cur_hp_len++;
    }
    hps[cur_pos] = cur_hp_start;      
    last_base = cur_base;
  }

  /* The last homopolymer back-fill hasn't been triggered
     yet by a new homopolymer, so we'll do it now...*/
  for( back_pos = cur_pos -1;
       back_pos >= cur_hp_start;
       back_pos-- ) {
    hpl[back_pos] = cur_hp_len;
  }
}

/* pop_s2c_in_a
   Args: (1) AlignmentP a - has a->seq2 and a->len2 set to 
   valid values
   Returns: void
   Populates the a->s2c array with code for quick lookup
   in submat
*/
void pop_s2c_in_a ( AlignmentP a ) {
  size_t i;
  int s_len;
  char b;
  s_len = a->len2;

  for ( i = 0; i < s_len; i++ ) {
    b = a->seq2[i];
    switch( b ) {
    case 'A' :
      a->s2c[i] = 0;
      break;
    case 'C' :
      a->s2c[i] = 1;
      break;
    case 'G' :
      a->s2c[i] = 2;
      break;
    case 'T' :
      a->s2c[i] = 3;
      break;
    default:
      a->s2c[i] = 4;
    }
  }
}
   

/* Input is a pointer to a valid Alignment. The value
   in a->len1 must be valid.
   Searches the last row (a->len1 -1) along all columns
   to find the best score that aligns all of the a->seq2
   Sets a->aec and a->aer to the correct values and sets
   a->best_score to the best score.
   Returns the best score */
int max_sg_score ( AlignmentP a ) {
  int row, col;
  int best_score = INT_MIN;
  row = a->len2 - 1;
  
  if ( row < 0 ) {
    /* Yeah, it happens if the sequence length is 0 */
    return best_score;
  }

  /* When finding the max_score, keep the earlier best score in
     case of a tie.
     This ensures that an alignment that is
     only against end-wrapped sequence will not be used if there
     is the same alignment earlier */
  for ( col = 0; col < a->len1; col++ ) {
    if ( a->m->mat[row][col].score > best_score ) {
      a->aec = col;
      a->aer = row;
      best_score = a->m->mat[row][col].score;
    }
  }
  a->best_score = best_score;
  return best_score;
}

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
		AlignmentP align) {
  /* Set up align with sequences */
  int col, row;
  int max_score = INT_MIN;
  align->seq1 = frag_seq->seq;

  align->len1 = strlen( align->seq1 );
  
  /* Populate the s1c and s2c arrays of submat indeces */
  pop_s1c_in_a( align );

  /* If a->hp then set up hpcl and hpcs for computing 
     special hp gap discount (hprl and hprs done already) */
  if ( align->hp ) {
    pop_hpl_and_hps( align->seq1, align->len1, 
		     align->hpcl, align->hpcs );
  }

  /* Align it! */
  dyn_prog( align );

  /* Good alignment? Check the last column and row of the
     dp matrix. 
     If the best score is >= TRIM_SCORE_CUT or if the number
     of bases trimmed x FLAT_MATCH = best_score, then set
     the trim point and trimmed flag in the frag_seq
  */
  col = align->len1 - 1;
  for( row = 0; row < align->len2; row++ ) {
    if ( align->m->mat[row][col].score > max_score ) {
      align->aec = col;
      align->aer = row;
      max_score = align->m->mat[row][col].score;
    }
  }

  /* Follow trace back to beginning of alignment */
  find_align_begin( align );

  /* Good enough? */
  if ( (max_score >= TRIM_SCORE_CUT) ||
       (max_score >= ((align->aer - align->abr + 1) 
		      * FLAT_MATCH)) ) {
    frag_seq->trimmed = 1;
    frag_seq->trim_point = align->abc - 1;
  }
  else {
    frag_seq->trimmed = 0;
  }
}

/* Takes pointers to two PWAlnFrag's
   The first one (front_pwaln) is populated by an alignment that
   crosses the wrap_point.
   Moves all of the alignment that is behind the wrap point into
   the back_pwaln and copies over all the other info
   Sets the correct segment flag for front and back */
void split_pwaln (PWAlnFragP front_pwaln, PWAlnFragP back_pwaln,
		  int wrap_point ) {
  int aln_pos; // current position in alignment
  int ref_pos; // current position in the reference sequence
  int frag_pos;// current position in the fragment sequence
  int aln_len; // length of reference alignment string

  aln_len = strlen( front_pwaln->ref_seq );
  ref_pos = front_pwaln->start;
  frag_pos = 0;
  aln_pos = 0;

  // Change the names of both fragments to avoid name conflicts with other tools
  int i = 0;
  while ((i <= (MAX_ID_LEN - 2)) && (front_pwaln->frag_id[i] != '\0'))
    i++; // Search for the end of fragment_id

  front_pwaln->frag_id[i]   = '_';
  front_pwaln->frag_id[i+1] = 'f';
  front_pwaln->frag_id[i+2] = '\0';

  strcpy( back_pwaln->frag_id, front_pwaln->frag_id );
  back_pwaln->frag_id[i+1] = 'b';    

  while ( ref_pos < wrap_point ) {
    /* Valid base at this position in ref_seq */
    if ( front_pwaln->ref_seq[aln_pos] != '-' ) {
      ref_pos++;
    }
    if ( front_pwaln->frag_seq[aln_pos] != '-' ) {
      frag_pos++;
    }
    aln_pos++;
  }

  /* OK, now aln_pos is at the first position that should go into
     the back_pwaln */
  strcpy( back_pwaln->ref_seq,  &front_pwaln->ref_seq[aln_pos] );
  strcpy( back_pwaln->frag_seq, &front_pwaln->frag_seq[aln_pos] );

  front_pwaln->ref_seq[aln_pos] = '\0';
  front_pwaln->frag_seq[aln_pos] = '\0';

  back_pwaln->start = 0;
  back_pwaln->end = front_pwaln->end;

  front_pwaln->end = wrap_point - 1;

  back_pwaln->segment = 'b';
  front_pwaln->segment = 'f';

  back_pwaln->offset = frag_pos;

  strcpy( back_pwaln->ref_id, front_pwaln->ref_id );
  strcpy( back_pwaln->ref_desc, front_pwaln->ref_desc );
//  strcpy( back_pwaln->frag_id, front_pwaln->frag_id );
  strcpy( back_pwaln->frag_desc, front_pwaln->frag_desc );

  back_pwaln->revcom     = front_pwaln->revcom;
  back_pwaln->trimmed    = front_pwaln->trimmed;
  back_pwaln->score      = front_pwaln->score;
  back_pwaln->num_inputs = front_pwaln->num_inputs;
}

int populate_pwaln_to_begin( AlignmentP a, PWAlnFragP pwaln ) {
  int row, col, next_row, next_col, ras_i, fas_i;
  char ras[ (INIT_ALN_SEQ_LEN * 2) + 1 ]; // temp place for constructing reference
  // alignment string
  char fas[ (INIT_ALN_SEQ_LEN * 2) + 1 ]; // temp place for constructing fragment
  //alignment string;

  ras[(INIT_ALN_SEQ_LEN*2)] = '\0';
  ras_i = (INIT_ALN_SEQ_LEN*2)-1;
  fas[(INIT_ALN_SEQ_LEN*2)] = '\0';
  fas_i = (INIT_ALN_SEQ_LEN*2)-1;

  row = a->aer;
  col = a->aec;

  while( (a->m->mat[row][col].trace != col) &&
	 (a->m->mat[row][col].trace != -row) ) {
    ras[ras_i--] = a->seq1[col];
    fas[fas_i--] = a->seq2[row];
    
    if ( a->m->mat[row][col].trace == 0 ) {
      row--;
      col--;
    }
      
    else {
	if ( a->m->mat[row][col].trace < 0 ) {
	  /* Negative number means gap up rows. So, cover sequence
	     of the fragment, but gaps in the reference */
	  next_row = -(a->m->mat[row][col].trace);
	  row--;
	  col--;
	  while( row > next_row ) {
	    fas[fas_i--] = a->seq2[row--];
	    ras[ras_i--] = '-';
	  }
	}
	else {
	  /* Positive number means gap back columns. Cover sequence
	     of the reference, but gaps in the fragment */
	  next_col = (a->m->mat[row][col].trace);
	  row--;
	  col--;
	  while( col > next_col ) {
	    fas[fas_i--] = '-';
	    ras[ras_i--] = a->seq1[col--];
	  }
	}
      }
    }

    ras[ras_i] = a->seq1[col];
    fas[fas_i] = a->seq2[row];
      
    strcpy( pwaln->ref_seq, &ras[ras_i] );
    strcpy( pwaln->frag_seq, &fas[fas_i] );
    return 1;
}


int sg_align ( MapAlignmentP maln, FragSeqP fs, FSDB fsdb, 
	       AlignmentP fw_a, AlignmentP rc_a, 
	       PWAlnFragP front_pwaln, 
	       PWAlnFragP back_pwaln) {
  int max_fw_score = INT_MIN;
  int max_rc_score = INT_MIN;
  RefSeqP rs;
  rs = maln->ref;
  AlignmentP best_a;

  fw_a->seq2 = fs->seq;

  rc_a->seq2 = fs->seq;
  
  if ( fs->trimmed ) {
    fw_a->len2 = fs->trim_point + 1;
    rc_a->len2 = fs->trim_point + 1;
  }
  else {
    fw_a->len2 = fs->seq_len;
    rc_a->len2 = fs->seq_len;
  }

  /* Populate the fw_a->s2c and rc_a->s2c codes */
  pop_s2c_in_a( fw_a );
  pop_s2c_in_a( rc_a );

  if ( fw_a->hp ) {
    pop_hpl_and_hps( fw_a->seq2, fw_a->len2,
		     fw_a->hprl, fw_a->hprs );
    pop_hpl_and_hps( rc_a->seq2, rc_a->len2,
		     rc_a->hprl, rc_a->hprs );
  }

  /* Set for a semiglobal alignment */
  fw_a->sg5 = 1;
  fw_a->sg3 = 1;
  rc_a->sg5 = 1;
  rc_a->sg3 = 1;

  /* Align it! */
  dyn_prog( fw_a );
  dyn_prog( rc_a );

  /* Find the best score */
  max_fw_score = max_sg_score( fw_a );
  max_rc_score = max_sg_score( rc_a );

  /* Which alignment has better score? */
  if ( max_fw_score > max_rc_score ) {
    best_a = fw_a;
  }
  else {
    best_a = rc_a;
  }

  find_align_begin( best_a );
 
  /* Load up front_pwaln */
  strcpy( front_pwaln->ref_id, rs->id );
  strcpy( front_pwaln->ref_desc, rs->desc );
  
  strcpy( front_pwaln->frag_id, fs->id );
  strcpy( front_pwaln->frag_desc, fs->desc );

  /* First, put all of alignment in front_pwaln */
  populate_pwaln_to_begin( best_a, front_pwaln );
      
  front_pwaln->start = best_a->abc;
  front_pwaln->end   = best_a->aec;
  front_pwaln->trimmed = fs->trimmed;
  front_pwaln->segment = 'a';
  front_pwaln->score = best_a->best_score;
  fs->score = best_a->best_score;

  /* Was this the rc alignment? */
  if ( best_a->rc ) {
    /* reverse complement the sequences */
    revcom_PWAF( front_pwaln );
    front_pwaln->revcom = 1;
    fs->rc = 1;
  }
  else {
    front_pwaln->revcom = 0;
    fs->rc = 0;
  }

  if ( best_a->rc ) {
    /* Adjust start and end coordinates if it's an rc alignment */
    front_pwaln->start = c2rcc( best_a->aec,
				rs->seq_len );
    front_pwaln->end   = c2rcc( best_a->abc,
				rs->seq_len );
    fs->as = c2rcc( best_a->aec, rs->seq_len );
    fs->ae = c2rcc( best_a->abc, rs->seq_len );
  }
  else {
    fs->as = best_a->abc;
    fs->ae = best_a->aec;
  }
  if ( fs->as > fs->ae ) {
    /* The fs->ae *should* be longer than the rs->seq_len
       for ease next time aligning */
    fs->ae = rs->seq_len + fs->as;
  }

  if ( front_pwaln->end > rs->seq_len ) {
    /* This alignment wraps around - adjust the end to
       demonstrate this for split_maln check */
    front_pwaln->end = front_pwaln->end - rs->seq_len;
  }

  /* Quit now if score is not good enough and distant_ref is not
     true */
  if ( (fs->score >= FIRST_ROUND_SCORE_CUTOFF) ||
       maln->distant_ref ) {
    
    /* Now, split up front_pwaln if this wrapped around the
       wrap point and merge them into the maln */
    if ( front_pwaln->start > front_pwaln->end ) {
      split_pwaln( front_pwaln, back_pwaln, rs->seq_len );
      if ( merge_pwaln_into_maln( front_pwaln, maln ) == 0 ) {
	return 0;
      }
      /* Point this fs->front_asp to the 
	 newly created AlnSeqP in maln */
      fs->front_asp = maln->AlnSeqArray[maln->num_aln_seqs - 1];

      if ( merge_pwaln_into_maln( back_pwaln, maln ) == 0 ) {
	return 0;
      }
      /* Point this fs->back_asp to the 
	 newly created AlnSeqP in maln */
      fs->back_asp = maln->AlnSeqArray[maln->num_aln_seqs - 1];
    }
    else {
      if ( merge_pwaln_into_maln( front_pwaln, maln ) == 0 ) {
	return 0;
      }
      /* Point this fs->front_asp to the 
	 newly created AlnSeqP in maln */
      fs->front_asp = maln->AlnSeqArray[maln->num_aln_seqs - 1];
      fs->back_asp = NULL;
    }

    /* Everyone is born unique until its discovered that they're not */
    fs->unique_best = 1;
    /* Every sequence has its own soul */
    fs->num_inputs  = 1;

    /* Did we see an alignment good enough for learning 
       what strand this is on, i.e., a positive-scoring 
       alignment? */
    if ( fs->score > FIRST_ROUND_SCORE_CUTOFF ) {
      fs->strand_known = 1;
    }
    else {
      fs->strand_known = 0;
    }

    if ( add_virgin_fs2fsdb( fs, fsdb ) == 0 ) {
      return 0;
    }
  }
  return 1;
}

