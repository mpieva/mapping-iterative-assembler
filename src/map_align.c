/* $Id: map_align.c 1111 2008-03-24 18:11:14Z green $ */
#include "map_align.h"


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

inline short int base2inx(const char base) {
	switch (base) {
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	default:
		return 4;
	}
}

int idCmp(const void* id1_, const void* id2_) {
	char** id1p = (char**) id1_;
	char** id2p = (char**) id2_;
	char* id1 = *id1p;
	char* id2 = *id2p;
	return strcmp(id1, id2);
}


/* Takes an AlnSeqP and a beginning and end coordinate of a region.
 All coordinates are 0-based
 Returns true is this AlnSeq overlaps the region at all, false
 if it does not */
inline int alnseq_ol_reg(AlnSeqP as, const int rs, const int re) {
    return ((as->start <= re) && (as->end >= rs));
}

/* This IDsList */
IDsListP init_ids_list(void) {
	IDsListP ids;
	char** ids_array;
	char* first_id;
	int i;

	// allocate the IDsList
	ids = (IDsListP)save_malloc(sizeof(IDsList));

	ids_array = (char**)save_malloc(INIT_NUM_IDS * sizeof( char* ));
	first_id = (char*)save_malloc(INIT_NUM_IDS * MAX_ID_LEN * sizeof(char));

	for (i = 0; i < INIT_NUM_IDS; i++) {
		ids_array[i] = &first_id[i*MAX_ID_LEN];
	}

	ids->num_ids = 0;
	ids->sorted = 0;
	ids->ids = ids_array;
	return ids;
}

void add_id(char* new_id, IDsListP used_ids_list) {
	/* No array out of bounds checking because we have a
	 hard limit for how many IDs we can ever see */
	strcpy(used_ids_list->ids[used_ids_list->num_ids], new_id);
	used_ids_list->num_ids++;
	qsort(used_ids_list->ids[0], used_ids_list->num_ids, (MAX_ID_LEN * sizeof(char)), idCmp);
	used_ids_list->sorted = 1;
}

void grow_ids_list(IDsListP ids) {
	int new_size, i, k;
	char** ids_array;
	char* first_id;
	new_size = (ids->size) * 2;

	ids_array = (char**)save_malloc(new_size * sizeof(char*));
	first_id = (char*)save_malloc(ids->size * MAX_ID_LEN * sizeof(char));

	/* Point first half of new ids_array to old half of pointers */
	for (i = 0; i < ids->size; i++) {
		ids_array[i] = ids->ids[i];
	}
	k = 0;
	/* Point secod half of new ids_array to new half of pointers */
	for (i = ids->size; i < new_size; i++) {
		ids_array[i] = &first_id[(k++ * MAX_ID_LEN)];
	}

	/* Free old ids */
	free(ids->ids);

	ids->ids = ids_array;
	ids->size = new_size;
}

int allowed_alignment(int ids_rest, IDsListP rest_ids_list, int no_dups,
		IDsListP used_ids_list, PWAlnFragP pwaln, double score_int,
		double score_slo) {
	char* found_id;
	int pass = 1; // Passes until we see a reason to exclude it

	/* Check the score */
	if (score_slo > -1) {
		if (pwaln->score < ((pwaln->end - pwaln->start + 1) * score_slo
				+ score_int)) {
			pass = 0;
		}
	}

	/* Score must be positive, no matter what! */
	if (pwaln->score <= 0) {
		pass = 0;
	}

	/* If we're restricted to a specified list of IDs, this one must
	 be on that list */
	if (ids_rest) {
		found_id = bsearch(pwaln->frag_id, rest_ids_list->ids[0],
				rest_ids_list->num_ids, (MAX_ID_LEN * sizeof(char)), idCmp);
		if ( !found_id) {
			pass = 0;
		}
	}

	/* If we won't allow duplicate IDs, this ID must not be a 
	 duplicate */
	if (no_dups) {
		found_id = bsearch(strcat(pwaln->frag_id, &pwaln->segment),
				used_ids_list->ids[0], rest_ids_list->num_ids,
				(MAX_ID_LEN * sizeof(char)), idCmp);
		if ( !(found_id == NULL)) {
			fprintf( stderr, "Already seen alignment for seq with ID %s\n",
			pwaln->frag_id);
			pass = 0;
		}
	}

	return pass;

}

int find_phred_qscore( BaseCountsP bcs ) {
  size_t i;
  int best_score;
  int not_best_scores[3];
  double p_best_score;
  double p_nbs[3];
  double p_correct;
  /* Is A the best-scoring base? */
  if ( (bcs->scoreA >= bcs->scoreC) &&
       (bcs->scoreA >= bcs->scoreG) &&
       (bcs->scoreA >= bcs->scoreT) ) {
    best_score = bcs->scoreA;
    not_best_scores[0] = bcs->scoreC;
    not_best_scores[1] = bcs->scoreG;
    not_best_scores[2] = bcs->scoreT;
  }

  else {
    /* Is C the best-scoring base? */
    if ( (bcs->scoreC >= bcs->scoreG) &&
	 (bcs->scoreC >= bcs->scoreT) ) {
      best_score = bcs->scoreC;
      not_best_scores[0] = bcs->scoreA;
      not_best_scores[1] = bcs->scoreG;
      not_best_scores[2] = bcs->scoreT;
    }
    else {
      /* Is G the best-scoring base? */
      if ( bcs->scoreG >= bcs->scoreT ) {
	best_score = bcs->scoreG;
	not_best_scores[0] = bcs->scoreA;
	not_best_scores[1] = bcs->scoreC;
	not_best_scores[2] = bcs->scoreT;
      }
      else {
	/* T must be the best-scoring base */
	best_score = bcs->scoreT;
	not_best_scores[0] = bcs->scoreA;
	not_best_scores[1] = bcs->scoreC;
	not_best_scores[2] = bcs->scoreG;
      }
    }
  }

  p_best_score = pow(2,((double)best_score/100));
  for( i = 0; i < 3; i++ ) {
    p_nbs[i] = pow(2,((double)not_best_scores[i]/100));
  }
  p_correct = p_best_score / ( p_nbs[0] + p_nbs[1] + p_nbs[2] );
  /* Check for overflow */
  if ( p_correct >= DBL_MAX ) {
    p_correct = DBL_MAX;
  }
  return 10 * log10(p_correct);
}

void show_single_pos(int ref_pos, char ref_base, 
		     char cons_base, BaseCountsP bcs) {
  printf("%d %c %c %d %d %d %d %d %d %d %d %d %d %d %0.3f\n", 
	 ref_pos, 
	 ref_base, 
	 cons_base,
	 bcs->cov, 
	 bcs->As, 
	 bcs->Cs, 
	 bcs->Gs, 
	 bcs->Ts, 
	 bcs->gaps,
	 bcs->scoreA,
	 bcs->scoreC,
	 bcs->scoreG,
	 bcs->scoreT,
	 find_phred_qscore( bcs ),
	 bcs->frac_agree
	 );
}

void add_base(char b, BaseCountsP bcs, PSSMP psm, int pssm_code) {
  short int b_inx;
  int depth;
  switch (b) {
  case 'A':
    bcs->As++;
    break;
  case 'C':
    bcs->Cs++;
    break;
  case 'G':
    bcs->Gs++;
    break;
  case 'T':
    bcs->Ts++;
    break;
  case '-':
    bcs->gaps++;
    break;
  }
  bcs->cov++;
  
  if (b == '-') {
    return;
  }
  
  b_inx = base2inx(b);
  depth = (pssm_code - 'A');
  
  bcs->scoreA += psm->sm[depth][0][b_inx];
  bcs->scoreC += psm->sm[depth][1][b_inx];
  bcs->scoreG += psm->sm[depth][2][b_inx];
  bcs->scoreT += psm->sm[depth][3][b_inx];
  
}

void reset_base_counts(BaseCountsP bc) {
  bc->As = 0;
  bc->Cs = 0;
  bc->Gs = 0;
  bc->Ts = 0;
  bc->scoreA = 0;
  bc->scoreC = 0;
  bc->scoreG = 0;
  bc->scoreT = 0;
  bc->gaps = 0;
  bc->cov = 0;
}

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
   Sets bcs->frac_agree to be the fraction of bases that are in 
   agreement with the highest-scoring base
*/
char find_consensus(BaseCountsP bcs, int cons_code) {
  int max_score;
  int top2scores[2];
  char max_base;

  /* Is this a zero coverage position that we really know 
     nothing about? */
  if (bcs->cov == 0) {
    bcs->frac_agree = 0.0;
    return 'N';
  }
  
  /* Is this a position whose consensus is a gap, i.e., no
     base at this position? */
  if ( ( (double)bcs->gaps / (double)bcs->cov ) >= (double)(PERC4GAP/100.0) ) {
    bcs->frac_agree = ((double)bcs->gaps / (double)bcs->cov);
    return '-';
  }
  

  /* If none of the above, find the highest and 2nd highest
     scoring base */
  top2scores[0] = bcs->scoreA;
  top2scores[1] = INT_MIN;
  max_score = bcs->scoreA;
  max_base = 'A';
  bcs->frac_agree = ((double)bcs->As / (double)bcs->cov);
  
  if (bcs->scoreC >= top2scores[0]) {
    max_score = bcs->scoreC;
    top2scores[1] = top2scores[0];
    top2scores[0] = bcs->scoreC;
    max_base = 'C';
    bcs->frac_agree = ((double)bcs->Cs / (double)bcs->cov);
  }
  else {
    top2scores[1] = bcs->scoreC;
  }
  
  /* Compare G to the best */
  if (bcs->scoreG >= top2scores[0]) {
    top2scores[1] = top2scores[0];
    top2scores[0] = bcs->scoreG;
    bcs->frac_agree = ((double)bcs->Gs / (double)bcs->cov);
    max_base = 'G';
  }
  else {
    if (bcs->scoreG >= top2scores[1]) {
      top2scores[1] = bcs->scoreG;
    }
  }

  /* Compare T to the best */
  if (bcs->scoreT >= top2scores[0]) {
    top2scores[1] = top2scores[0];
    top2scores[0] = bcs->scoreT;
    max_base = 'T';
    bcs->frac_agree = ((double)bcs->Ts / (double)bcs->cov);
  }
  else {
    if (bcs->scoreT >= top2scores[1]) {
      top2scores[1] = bcs->scoreT;
    }
  }
  
  if ( cons_code == 1 ) {
    /* If the best scoring base was better than the
       MIN_SCORE_CON, that is the consensus. 
       If none were, consensus is N */
    if (top2scores[0] >= MIN_SCORE_CONS) {
      return max_base;
    } 
    else {
      return 'N';
    }
  }

  if ( cons_code == 2 ) {
    if ((top2scores[0] >= 0) || 
	(top2scores[0] - MIN_SC_DIFF_CONS) > top2scores[1] ) {
      return max_base;
    }
    else {
      return 'N';
    }
  }

  /* Default => act like cons_code == 1 */

  else {
    if (top2scores[0] >= MIN_SCORE_CONS) {
      return max_base;
    } 
    else {
      return 'N';
    }
  }
}

int alnSeqCmp(const void* as1_, const void* as2_) {
	const AlnSeqP *as1 = (const AlnSeqP*) as1_ ;
	const AlnSeqP *as2 = (const AlnSeqP*) as2_ ;
	if ( ((*as1)->start) < ((*as2)->start)) {
		return -1;
	}
	if ( ((*as1)->start) > ((*as2)->start)) {
		return 1;
	}
	if ( ((*as1)->start) == ((*as2)->start)) {
		if ( ((*as1)->end) < ((*as2)->end)) {
			return -1;
		}
		if ( ((*as1)->end) > ((*as2)->end)) {
			return 1;
		}
		if ( ((*as1)->end) == ((*as2)->end)) {
			return 0;
		}
	}
	return 0;
}

/* Computes the reverse complement of a single base.  It must support
 * IUPAC ambiguity codes, small case letters, and gap symbols. */
char revcom_char(const char base) {
    char tbl[] = "TVGH\0\0CD\0\0M\0KN\0\0\0YSAABWXR\0" ;
    char rc = 0 ;

    if( base == '-' ) return '-' ;

    if( 'A' <= base && base <= 'Z' )
        rc = tbl[ base - 'A' ] ;
    else if( 'a' <= base && base <= 'z' )
        rc = tbl[ base - 'a' ] + 32 ;

    if( rc ) return rc ;
    fprintf( stderr, "Do not know how to revcom \"%c\"\n", base);
    return 'N';
}

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
		int out_format) {
	int i, j, ins_len, this_frag_ins_len;
	char* ins_seq;
	AlnSeqP aln_seq;
	BaseCountsP* bcs_array;
	BaseCountsP first_bcs;
	PSSMP psm;

	ins_len = maln->ref->gaps[pos];

	bcs_array = (BaseCountsP*)save_malloc(ins_len * sizeof(BaseCountsP));
	first_bcs = (BaseCountsP)save_malloc(ins_len * sizeof(BaseCounts));

	for (i = 0; i < ins_len; i++) {
		bcs_array[i] = &first_bcs[i];
		reset_base_counts(bcs_array[i]);
	}

	for (i = 0; i < maln->num_aln_seqs; i++) {
		aln_seq = maln->AlnSeqArray[i];
		/* Does this aligned fragment cover this position? */
		if ( (aln_seq->start < pos) && // It does not cover this position
				//if it starts exactly here because the gap is, by convention,
				//just upstream of this position
				(aln_seq->end >= pos)) {
			if (aln_seq->revcom) {
				psm = maln->rpsm;
			} else {
				psm = maln->fpsm;
			}
			/* Does it have some actual inserted sequence? */
			ins_seq = aln_seq->ins[pos - aln_seq->start];
			if (ins_seq == NULL) {
				for (j = 0; j < ins_len; j++) {
					add_base( '-', bcs_array[j], psm,
							aln_seq->smp[pos - aln_seq->start]);
				}
			} else {
				this_frag_ins_len = strlen(ins_seq);
				for (j = 0; j < ins_len; j++) {
					if (j < this_frag_ins_len) {
						add_base(ins_seq[j], bcs_array[j], psm,
								aln_seq->smp[pos - aln_seq->start]);
					} else {
						add_base( '-', bcs_array[j], psm,
								aln_seq->smp[pos - aln_seq->start]);
					}
				}
			}
		}
	}

	for (j = 0; j < ins_len; j++) {
		ins_cons[j] = find_consensus(bcs_array[j], maln->cons_code);
		cons_cov[j] = bcs_array[j]->cov;
		if ( (out_format == 4) && !(ins_cons[j] == '-')) {
			show_single_pos(pos, '-', ins_cons[j], bcs_array[j]);
		}
		if (out_format == 41) {
			show_single_pos(pos, '-', ins_cons[j], bcs_array[j]);
		}
	}

	free(first_bcs);
	free(bcs_array);
}

void revcom_PWAF(PWAlnFragP pwaln) {
  char tmp_ref, tmp_frag;
  int len, i;
  len = strlen(pwaln->ref_seq);

  for (i = 0; i < len/2; i++) {
    tmp_ref = pwaln->ref_seq[i];
    pwaln->ref_seq[i] = revcom_char(pwaln->ref_seq[len-(i+1)]);
    pwaln->ref_seq[len-(i+1)] = revcom_char(tmp_ref);
    
    tmp_frag = pwaln->frag_seq[i];
    pwaln->frag_seq[i] = revcom_char(pwaln->frag_seq[len-(i+1)]);
    pwaln->frag_seq[len-(i+1)] = revcom_char(tmp_frag);
  }
  
  /* If sequence length is even, we're done, otherwise there is
     the base right in the center to revcom */
  if (len%2 == 1) {
    pwaln->ref_seq[i] = revcom_char(pwaln->ref_seq[len-(i+1)]);
    pwaln->frag_seq[i] = revcom_char(pwaln->frag_seq[len-(i+1)]);
  }
  pwaln->revcom = 1;

}


/* For a given region, defined by reg_start and reg_end, show
 the refence sequence, the consensus sequence, 
 and the sequence of all the fragments that overlap this
 region at all.
 */
void print_region( MapAlignmentP maln, int reg_start, int reg_end,
		   int out_format, int in_color ) {
  int i, ref_pos, ref_gaps, j, cons_pos, ins_len;
  int num_gaps = 0;
  int ins_seq_len;
  int read_out_pos;
  char* consensus;
  char* aln_ref;
  char* read_reg;
  char* ins_cons;
  char* read_str;
  char* read_id;
  char* ins_seq;
  int* ins_cov;
  BaseCountsP bcs;
  AlnSeqP aln_seq;
  PSSMP psm;
  
  /* Make sure region doesn't go off edge */
  if (reg_start < 1) {
    reg_start = 1;
  }
  if (reg_end > maln->ref->seq_len) {
    reg_end = maln->ref->seq_len;
  }
  
  bcs = (BaseCountsP)save_malloc(sizeof(BaseCounts));
  reset_base_counts(bcs);
  
  /* Find how many gaps are in this region */
  for (i = reg_start-1; i <= reg_end; i++) {
    num_gaps += maln->ref->gaps[i];
  }
  
  /* Make char arrays long enough for the sequence plus
     gaps for the reference, the consensus, and a single 
     read. These will be populated and output by the rest
     of this function.
  */
  consensus = (char*)save_malloc((num_gaps + (reg_end-reg_start+1) + 10)
				 * sizeof(char));
  aln_ref = (char*)save_malloc((num_gaps + (reg_end-reg_start+1) + 10)
			       * sizeof(char));
  read_reg = (char*)save_malloc((num_gaps + (reg_end-reg_start+1) + 10)
				* sizeof(char));
	
  /* Make char and int array for insert consensus and
     insert coverage to be used whenever needed */
  ins_cons = (char*)save_malloc(MAX_INS_LEN * sizeof(char));
  ins_cov = (int* )save_malloc(MAX_INS_LEN * sizeof(int));
  
  cons_pos = 0;
  for (ref_pos = reg_start - 1; ref_pos < reg_end; ref_pos++) {
    ref_gaps = maln->ref->gaps[ref_pos];
    /* Add these gaps to the reference aligned string and the inserted
       sequence to the consensus[] */
    if (ref_gaps > 0) {
      find_ins_cons(maln, ref_pos, ins_cons, ins_cov, out_format);
      for (j = 0; j < ref_gaps; j++) {
	aln_ref[cons_pos] = '-';
	consensus[cons_pos] = ins_cons[j];
	cons_pos++;
      }
    }
    /* Re-zero all the base counts */
    reset_base_counts(bcs);
    
    /* Find all the aligned fragments that include this
       position and make a consensus from it */
    for (j = 0; j < maln->num_aln_seqs; j++) {
      aln_seq = maln->AlnSeqArray[j];
      /* Does this aligned fragment cover this position? */
      if ( (aln_seq->start <= ref_pos) && // checked
	   (aln_seq->end >= ref_pos)) {
	if (aln_seq->revcom) {
	  psm = maln->rpsm;
	} else {
	  psm = maln->fpsm;
	}
	add_base(aln_seq->seq[ref_pos - aln_seq->start], bcs, psm,
		 aln_seq->smp[ref_pos - aln_seq->start]);
      }
    }
    
    consensus[cons_pos] = find_consensus(bcs, maln->cons_code);
    aln_ref[cons_pos] = maln->ref->seq[ref_pos];
    cons_pos++;
  }
  
  consensus[cons_pos] = '\0';
  aln_ref[cons_pos] = '\0';
  
  /* Now print the reference and the consensus */
  if (out_format == 61) {
    fasta_aln_print(aln_ref, maln->ref->id);
    fasta_aln_print(consensus, "Consensus");
  } else {
    if (in_color) {
      printf("%-20.20s ", maln->ref->id);
      color_print(aln_ref);
      printf("%-20.20s ", "Consensus");
      color_print(consensus);
    } else
      printf("%-20.20s %s\n%-20s %s\n", maln->ref->id, aln_ref,
	     "Consensus", consensus);
  }
  
  /* 
     Alloc memories for the string to hold each read (plus .'s outside)
     and alloc memories for the special id which is the regular ID
     plus the code for whether it's truncated, reverse complemented,
     and the number of input sequence
  */
  read_str = (char*)save_malloc(strlen(aln_ref) * sizeof(char) + 1);
  read_id  = (char*)save_malloc((MAX_ID_LEN + 4) * sizeof(char) + 1);
  /* Find every sequence that overlaps this region and print
     the overlapping segment */
  for (j = 0; j < maln->num_aln_seqs; j++) {
    aln_seq = maln->AlnSeqArray[j];
    if (alnseq_ol_reg(aln_seq, (reg_start-1), (reg_end-1)) ) {
      read_out_pos = 0;
      if (aln_seq->trimmed) {
	read_id[0] = 't';
      } else {
	read_id[0] = '_';
      }
      
      if (aln_seq->revcom) {
	read_id[1] = 'r';
      } else {
	read_id[1] = '_';
      }
      sprintf( &read_id[2], "%02d", aln_seq->num_inputs );
      read_id[4] = '\0';
           
      char * out_read_id = calloc(strlen(read_id) + strlen(aln_seq->id)+1, sizeof(char));
      strcpy(out_read_id, aln_seq->id);
      strcat(out_read_id, read_id);      
      out_read_id[strlen(read_id) + strlen(aln_seq->id)] = '\0';
          
      
      if (out_format == 6) {        
	printf("%-20.20s ", out_read_id);
      }
      for (ref_pos = reg_start - 1; ref_pos < reg_end; ref_pos++) {
	ref_gaps = maln->ref->gaps[ref_pos];
	/* Check to make sure that this fragment has started and
	   not ended by this ref_pos */
	if ( (aln_seq->start <= ref_pos) && // checked
	     (aln_seq->end >= ref_pos)) {
	  if (ref_gaps > 0) {
	    if (aln_seq->ins[ref_pos - aln_seq->start] == NULL) {
	      ins_len = 0;
	    } else {
	      ins_len
		= strlen(aln_seq->ins[ref_pos - aln_seq->start]);
	    }
	    if (aln_seq->start == ref_pos) {
	      // Exactly at the beginning of this frag
	      for (i = 0; i < ref_gaps; i++) {
		read_str[read_out_pos++] = '.';
		//		printf( "." );
	      }
	    } else {
	      // Just a normal, interior gapped position
	      if (ins_len > 0) {
		ins_seq
		  = aln_seq->ins[ref_pos - aln_seq->start];
		ins_seq_len = strlen(ins_seq);
		for (i = 0; i < ins_seq_len; i++) {
		  read_str[read_out_pos++] = ins_seq[i];
		}
		//		printf( "%s", aln_seq->ins[ref_pos - aln_seq->start] );
	      }
	      for (i = 0; i < (ref_gaps - ins_len); i++) {
		read_str[read_out_pos++] = '-';
		//		printf( "-" );
	      }
	    }
	  }
	  read_str[read_out_pos++]
	    = aln_seq->seq[ref_pos - aln_seq->start];
	  //printf( "%c", aln_seq->seq[ref_pos - aln_seq->start] );
	} else {
	  // This fragment doesn't actually cover this base
	  for (i = 0; i < ref_gaps; i++) {
	    // print this . for all ref gaps
	    read_str[read_out_pos++] = '.';
	    // printf( "." );
	  }
	  read_str[read_out_pos++] = '.';
	  //printf( "." );
	}
      }
      read_str[read_out_pos] = '\0';
      if (out_format == 61) {
	fasta_aln_print(read_str, out_read_id);
      } else {
	
	if (in_color) {
	  color_print(read_str);
	} else
	  printf("%s\n", read_str);
      }
      
      free(out_read_id);
    }
  }
  free(bcs);
  free(consensus);
  free(aln_ref);
  free(read_reg);
  free(ins_cons);
  free(ins_cov);
  free(read_str);
  free(read_id);
}

void col_print_cons(char* consensus, char* aln_ref, int* cov, int* ref_poss,
		MapAlignmentP maln) {
	int len, i;
	char c;
	int* starts_f;
	int* starts_r;
	int* ends_f;
	int* ends_r;
	AlnSeqP as;

	len = strlen(consensus);

	starts_f = (int* )save_malloc(len * sizeof(int));
	starts_r = (int* )save_malloc(len * sizeof(int));
	ends_f = (int* )save_malloc(len * sizeof(int));
	ends_r = (int* )save_malloc(len * sizeof(int));

	/* Initialize everything to zero */
	for (i = 0; i < len; i++) {
		starts_f[i] = 0;
		starts_r[i] = 0;
		ends_f[i] = 0;
		ends_r[i] = 0;
	}

	/* Now, go through all the aligned fragments and update
	 the starts and ends arrays based on where each fragment...
	 starts and ends! */
	for (i = 0; i < maln->num_aln_seqs; i++) {
		as = maln->AlnSeqArray[i];
		if (as->revcom) {
			switch (as->segment) {
			case 'f':
				/* only the start is correct for front fragments */
				starts_r[as->start]++;
				break;
			case 'b':
				/* only the end is correct for back fragments */
				ends_r[as->end]++;
				break;
			default:
				starts_r[as->start]++;
				ends_r[as->end]++;
				break;
			}
		}

		/* Not reverse complement */
		else {
			switch (as->segment) {
			case 'f':
				/* only the start is correct for front fragments */
				starts_f[as->start]++;
				break;
			case 'b':
				/* only the end is correct for back fragments */
				ends_f[as->end]++;
				break;
			default:
				starts_f[as->start]++;
				ends_f[as->end]++;
				break;
			}
		}
	}

	printf("# Columns:\n");
	printf("# 1. Assembly consensus base\n");
	printf("# 2. Reference %s base\n", maln->ref->id);
	printf("# 3. Coverage (number of reads overlapping this position)\n");
	printf("# 4. Coordinate on reference sequence (1-based)\n");
	printf("# 5. Number of fragments on forward strand that start here\n");
	printf("# 6. Number of fragments on reverse strand that start here\n");
	printf("# 7. Number of fragments on forward strand that end here\n");
	printf("# 8. Number of fragments on reverse strand that end here\n");
	for (i = 0; i < len; i++) {
		if ( !((consensus[i] == '-') && (aln_ref[i] == '-') )) {
			if (consensus[i] == ' ') {
				c = 'X';
			} else {
				c = consensus[i];
			}
			printf("%c\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", c, aln_ref[i], cov[i],
					(ref_poss[i]+1), starts_f[ref_poss[i]],
					starts_r[ref_poss[i]], ends_f[ref_poss[i]],
					ends_r[ref_poss[i]]);
		}
	}
}



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
int merge_pwaln_into_maln(PWAlnFragP pwaln, MapAlignmentP maln) {
  int i, j, aln_len, ref_frag_len, ref_pos, gap_compare, mind_the_gap,
    seq_pos;
  char c, f;
  char* ins_seq;
  AlnSeqP asp;
  int this_ref_gaps[(2*INIT_ALN_SEQ_LEN) + 1];
  
  // Grow array of aligned sequences if necessary
  if (maln->num_aln_seqs >= maln->size) {
    if ( !(grow_alns_map_alignment(maln))) {
      return 0;
    }
  }
  
  // Get a pointer to this next AlnSeq
  asp = maln->AlnSeqArray[maln->num_aln_seqs];
  
  // Copy over all the details thusfar
  strcpy(asp->id, pwaln->frag_id);
  strcpy(asp->desc, pwaln->frag_desc);
  asp->score      = pwaln->score;
  asp->start      = pwaln->start;
  asp->end        = pwaln->end;
  asp->revcom     = pwaln->revcom;
  asp->trimmed    = pwaln->trimmed;
  asp->segment    = pwaln->segment;
  asp->num_inputs = pwaln->num_inputs;
  aln_len = strlen(pwaln->frag_seq);

  /* Copy the fragment aligned sequence string, gap characters
     and all, into asp->seq 
  */
  mind_the_gap = 0;
  j = 0;
  seq_pos = 0;
  this_ref_gaps[seq_pos] = 0;
  for (i = 0; i < aln_len; i++) {
    c = pwaln->ref_seq[i];
    f = pwaln->frag_seq[i];
    if (c == '-') {
      this_ref_gaps[seq_pos]++;
      if (mind_the_gap) {
	// Extending an already started gap
	ins_seq[j++] = pwaln->frag_seq[i];
      } else {
	// Starting a new gap
	ins_seq = (char*)save_malloc(MAX_INS_LEN * sizeof(char));
	j = 0;
	ins_seq[j++] = f;
      }
      mind_the_gap = 1;
    } 
    else {
      // Not a gap
      if (mind_the_gap) {
	// Just finished a gap, add \0 to inserted sequence
	ins_seq[j] = '\0';
	asp->ins[seq_pos] = ins_seq;
      } 
      else { // Not a gap here
	asp->ins[seq_pos] = NULL;
      }
      asp->seq[seq_pos++] = f;
      this_ref_gaps[seq_pos] = 0;
      mind_the_gap = 0;
    }
  }
  
  /* Add string terminator, just in case */
  asp->seq[seq_pos] = '\0';
  
  // Now, go through these ref seq gaps and see if they were already
  // known before
  ref_frag_len = asp->end - asp->start + 1;
  for (i = 0; i < ref_frag_len; i++) {
    ref_pos = asp->start + i;
    gap_compare = this_ref_gaps[i] - maln->ref->gaps[ref_pos];
    
    if (gap_compare > 0) {
      /* Longer gap in this fragment than known before so we must
	 make maln->ref->gaps[ref_pos] longer to accomodate it
      */
      maln->ref->gaps[ref_pos] += gap_compare;
    }
  }
  maln->num_aln_seqs++;
  return 1;
}



/* Takes the description from an Udo align aligned sequence
 and puts the start, end, strand, and score information in
 the correct field of the PWAlnFragP
 The desc (description) is a string like this, e.g.:
 "- 4199-4261 score=5441"
 */
int ses_from_align_desc(PWAlnFragP pwaln, int* strand) {
	//  char* desc, int* start,
	//		 int* end, int* strand ) {
	char strand_char;
	char score[16];
	int parts_read = 0;
	pwaln->segment = 'n'; // not applicable
	parts_read = sscanf(pwaln->ref_desc, "%c %d-%d score=%s %c", &strand_char,
			&pwaln->start, &pwaln->end, score, &pwaln->segment);

	if (parts_read < 4) {
		// failure
		return 0;
	}

	// Convert the 1-based coordinates to 0-based
	pwaln->start--;
	pwaln->end--;

	if (score[0] == '-') {
		pwaln->score = -1*(atoi(&score[1]));
	} else {
		pwaln->score = atoi(score);
	}

	if (strand_char == '+') {
		*strand = 1;
		return 1;
	}
	if (strand_char == '-') {
		*strand = -1;
		return 1;
	}

	return 0;
}

/* adapt_from_desc checks the frag_desc string of a PWAlnFrag,
 given a pointer to one (PWAlnFragP) and sets the trimmed
 field to true (1) if the phrase "adapter cut off" is there
 Returns true if everthing went fine, false otherwise
 */
int adapt_from_desc(PWAlnFragP af) {
	if (af->frag_desc == NULL) {
		return 0;
	}

	if (strstr(af->frag_desc, "adapter cut off") == NULL) {
		af->trimmed = 0;
	} else {
		af->trimmed = 1;
	}
	return 1;
}


/* Grow the space for a sequence (an array of char)
 to twice its current size
 Copy its current contents into the new sequence
 Free the now unused old memory
 */
char* grow_seq(char* seq, int size) {
  int i;
  char* new_seq;
  new_seq = (char*)save_malloc( 2 * size );
  for (i = 0; i < size; i++) {
    new_seq[i] = seq[i];
  }
  free(seq);
  return new_seq;
}


