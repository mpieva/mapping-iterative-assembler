#include "kmer.h"

/* kmer2inx
   Args: (1) a pointer to a character string;
             the kmer to find the corresponding index of;
	     might not be null-terminated
	 (2) length of the kmer
	 (3) pointer to size_t to put the index
   Returns: TRUE if the index was set, FALSE if it could not
            be set because of some non A,C,G,T character
   Uses the formula A=>00, C=>01, G=>11, T=>11 to make a
   bit string for the kmer. Any other character is not allowed
   and will cause an error
   The bit string is constructed by reading the kmer from left
   to right. This bit-string is then interpreted as a variable
   of type size_t and is appropriate as an array index
*/
int kmer2inx( const char* kmer,
	      const unsigned int kmer_len,
	      size_t* inx ) {
  size_t l_inx = 0;
  int i = 0;
  char curr_char;

  while( i < kmer_len ) {
    l_inx = l_inx << 2;
    curr_char = toupper(kmer[i]); // Upper case it in case it is not
    switch( curr_char ) {
    case 'A' :
      l_inx += 0;
      break;
    case 'C' :
      l_inx += 1;
      break;
    case 'G' :
      l_inx += 2;
      break;
    case 'T' :
      l_inx += 3;
      break;
    default :
      return 0; // not valid!
    }
    i++;
  }
  *inx = l_inx;
  return 1; // valid!
}


/* add_kmer
   Args: (1) KPL* kmer array
         (2) index position - must be valid
         (3) position of this kmer to add
   Returns: void
   Takes a newly discovered kmer position and adds it to the
   array. The given index specifies what the kmer is, but
   for this operation, we really do not care. We simply add
   the position to the positions field of this kmer. If this
   kmer has never been seen before, then we have to
   initialize it, too.
*/
void add_kmer( KPL* kpa, const size_t inx, const size_t i ) {

  /* Check if it is not already initialized */
  if ( kpa[inx] == NULL ) {
    /* Never seen this kmer before, so set it up */
    kpa[inx] = (KPL)save_malloc(sizeof(KmerPosList));
    kpa[inx]->num_pos = 0;
    kpa[inx]->sorted  = 0;
  }

  /* Check to make sure we're not over the maximum number
     of allowable positions for this kmer */
  if ( kpa[inx]->num_pos == MAX_KMER_POS ) {
    return;
  }

  /* No? then add it */
  kpa[inx]->positions[kpa[inx]->num_pos] = i;
  kpa[inx]->sorted = 0;
  kpa[inx]->num_pos++;

  return;
}
/* init_kpa
   Args: (1) length of kmers to use
   Returns: pointer to KPL; an array of pointers to KmerPosList
*/
KPL* init_kpa( const int kmer_len ) {
  KPL* kpa;
  unsigned int size = 1;
  if ( kmer_len > MAX_KMER_LEN ) {
    fprintf( stderr, "Cannot use kmer length greater than %d\n",
	     MAX_KMER_LEN );
    exit( 2 );
  }
  size = size << (2*kmer_len);
  kpa = (KPL*)calloc(size, sizeof(KPL));
  if ( kpa == NULL ) {
    fprintf( stderr,
	     "Not enough memories for kmers of length %d\n",
	     kmer_len );
    exit( 1 );
  }
  return kpa;
}


void grow_kmers ( KmersP k ) {
  int new_size, i, j;
  char** new_kmers;
  char* first_id;
  new_size = k->size * 2;
  new_kmers = (char**)save_malloc( new_size * sizeof(char*) );
  first_id  = (char*)save_malloc(k->size * (k->kmer_len + 1) * sizeof(char));

  /* Point first half of new_kmers to already existing pointers */
  for( i = 0; i < k->size; i++ ) {
    new_kmers[i] = k->kmers[i];
  }
  j = 0;
  /* Point the second half of new_kmers to new points */
  for( i = k->size; i < new_size; i++ ) {
    new_kmers[i] = &first_id[(k->kmer_len + 1) * j++];
  }

  /* Free old kmers */
  free( k->kmers );
  k->kmers = new_kmers;
  k->size = new_size;
}

/* all_upper
   Args: (1) Pointer to char array (seq)
         (2) int number of characters to check (len)
   Returns: int 1 => first len charaters in seq are all upper case
                0 => at least one of the characters is not upper case
*/
inline int all_upper( const char* seq, const int kmer_len ) {
  size_t i;
  for( i = 0; i < kmer_len; i++ ) {
    if ( islower( seq[i] ) ) {
      return 0;
    }
  }
  return 1;
}


/* populate_kpa
 */
int populate_kpa( KPL* kpa, const char* seq,
		  const size_t seq_len,
		  const int kmer_len,
		  const int soft_mask ) {
  size_t i, inx;
  for( i = 0; i <= (seq_len - kmer_len); i++ ) {
    /* Add this kmer if we're not check for softmasking or
       if we are and it passes the test */
    if ( !soft_mask || all_upper(&seq[i], kmer_len) ) {
      if ( kmer2inx( &seq[i], kmer_len, &inx ) ) {
	add_kmer( kpa, inx, i );
      }
    }
  }
  return 1;
}



/* pop_kmers
   Args: (1) RefSeqP ref - reference sequence with forward and reverse sequence
         (2) int kmer_filt_len - length of kmers
   Initializes a Kmers struct and populates it with all the kmers in the
   forward and reverse-complement sequence of the input RefSeq.
   Returns: pointer to KmersP
*/
KmersP pop_kmers( RefSeqP ref, int kmer_filt_len ) {
  int i, pos;
  KmersP k;
  char* first_kmer;
  char* curr_kmer;
  /* Allocate memory for the kmers, guessing how much will be required */
  k = (KmersP)save_malloc(sizeof( struct kmers ) );
  k->num_kmers = 0;
  k->kmer_len = kmer_filt_len;
  k->size = ref->seq_len; // just a guess
  k->kmers = (char**)save_malloc( k->size * sizeof( char* ) );
  first_kmer = (char*)save_malloc( k->size * (kmer_filt_len+1) * sizeof(char) );
  for( i = 0; i < k->size; i++ ) {
    k->kmers[i] = &first_kmer[i * (kmer_filt_len+1)];
  }

  curr_kmer = (char*)save_malloc( (kmer_filt_len + 1) * sizeof(char) );
  curr_kmer[kmer_filt_len] = '\0'; // Null-terminate for now and forever

  /* Now, cruise through the forward and reverse complement and
     load up the kmers */
  for( pos = 0; pos <= ref->wrap_seq_len - kmer_filt_len; pos++ ) {
    strncpy( curr_kmer, &ref->seq[pos], kmer_filt_len );
    if ( bsearch( &curr_kmer, k->kmers, k->num_kmers, sizeof(char*), idCmp )
	 == NULL ) {
      /* not found, add it */
      k->num_kmers++;
      if ( k->num_kmers >= k->size ) {
	grow_kmers( k );
      }
      strncpy( k->kmers[k->num_kmers - 1], curr_kmer, (kmer_filt_len+1) );
      qsort( k->kmers, k->num_kmers, sizeof(char*), idCmp );
    }
  }

  for( pos = 0; pos <= ref->wrap_seq_len - kmer_filt_len; pos++ ) {
    strncpy( curr_kmer, &ref->rcseq[pos], kmer_filt_len );
    if ( bsearch( &curr_kmer, k->kmers, k->num_kmers, sizeof(char*), idCmp )
	 == NULL ) {
      /* not found, add it */
      k->num_kmers++;
      if ( k->num_kmers >= k->size ) {
	grow_kmers( k );
      }
      strncpy( k->kmers[k->num_kmers], curr_kmer, (kmer_filt_len+1) );
      qsort( k->kmers, k->num_kmers, sizeof(char*), idCmp );
    }
  }

  k->sorted = 1;

  free( curr_kmer );
  return k;
}


/* Returns: TRUE (1) if we should align this sequence
            FALSE (0) if we should NOT align this sequence because
                      it shares no kmers with the reference
*/
int new_kmer_filter( FragSeqP fs,
		     KPL* fkpa,
		     KPL* rkpa,
		     int kmer_len,
		     AlignmentP fwa,
		     AlignmentP rca ) {
  size_t frag_len, frag_pos, inx, ref_len, ref_pos, i;
  int mask_min, mask_max; // Sometimes these become negative
  unsigned int num_f_kmers_found = 0;
  unsigned int num_r_kmers_found = 0;

  /* Check for no kmer filtering */
  if ( kmer_len < 0 ) {
    memset( fwa->align_mask, 1, fwa->len1 );
    memset( fwa->align_mask, 1, rca->len1 );
    return 1;
  }

  /* Reset the alignment masks */
  memset( fwa->align_mask, 0, fwa->len1 );
  memset( rca->align_mask, 0, rca->len1 );

  /* How long is this fragment? */
  if ( fs->trimmed ) {
    frag_len = (fs->trim_point + 1);
  }
  else {
    frag_len = fs->seq_len;
  }

  if ( frag_len < kmer_len ) {
    return 0;
  }

  /* Zip through all the kmers in this fragment sequence. If any
     are present in the forward or reverse kpa's, then we pass
     the filter, i.e., return 1 */
  for( frag_pos = 0; frag_pos <= (frag_len - kmer_len); frag_pos++ ) {
    if ( kmer2inx( &fs->seq[frag_pos], kmer_len, &inx ) ) {
      if ( fkpa[inx] != NULL ) {
	ref_len = fwa->len1;
	/* There are some kmers here. Add them to the total
	   count and update the align_mask */
	num_f_kmers_found += fkpa[inx]->num_pos;
	if ( num_f_kmers_found >= KMER_SATURATE ) {
	  memset( fwa->align_mask, 1, fwa->len1 );
	}

	for( i = 0; i < fkpa[inx]->num_pos; i++ ) {
	  /* Unmask the region surrounding this kmer */
	  ref_pos = fkpa[inx]->positions[i];
	  mask_min = ref_pos - frag_pos - ALIGN_MASK_BUFFER;
	  if ( mask_min < 0 ) {
	    mask_min = 0;
	  }
	  mask_max = ref_pos + (frag_len - frag_pos) + ALIGN_MASK_BUFFER;
	  if ( mask_max >= ref_len ) {
	    mask_max = (ref_len - 1);
	  }
	  memset( &fwa->align_mask[mask_min], 1, (mask_max-mask_min+1) );
	}
      }

      if ( rkpa[inx] != NULL ) {
	ref_len = rca->len1;
	/* There are some kmers here. Add them to the total
	   count and update the align_mask */
	num_r_kmers_found += rkpa[inx]->num_pos;
	if ( num_r_kmers_found >= KMER_SATURATE ) {
	  memset( rca->align_mask, 1, rca->len1 );
	}

	for( i = 0; i < rkpa[inx]->num_pos; i++ ) {
	  /* Unmask the region surrounding this kmer */
	  ref_pos = rkpa[inx]->positions[i];
	  mask_min = ref_pos - frag_pos - ALIGN_MASK_BUFFER;
	  if ( mask_min < 0 ) {
	    mask_min = 0;
	  }

	  mask_max = ref_pos + frag_len - frag_pos - 1 + ALIGN_MASK_BUFFER;
	  if ( mask_max >= ref_len ) {
	    mask_max = (ref_len - 1);
	  }
	  memset( &rca->align_mask[mask_min], 1, (mask_max-mask_min+1) );
	}
      }
    }
  }

  /* Return 0 if no kmers found; TRUE (not 0) if some kmers found */
  return (num_f_kmers_found + num_r_kmers_found);
}

int kmer_filter( int kmer_filt_len, FragSeqP fs, KmersP k ) {
  int len, pos;
  char* test_kmer;

  /* First, check if the user wants any kmer filtering
     Special value -1 means none and it doesn't make sense
     to filter for kmers <= 0 */
  if ( kmer_filt_len < 0 ) {
    return 1;
  }

  test_kmer = (char*)save_malloc( (k->kmer_len + 1) * sizeof(char) );
  test_kmer[k->kmer_len] = '\0';

  if ( fs->trimmed ) {
    len = fs->trim_point;
  }
  else {
    len = fs->seq_len - 1;
  }

  for( pos = 0; pos <= len - k->kmer_len; pos++ ) {
    strncpy( test_kmer, &fs->seq[pos], k->kmer_len );
    if ( bsearch( &test_kmer, k->kmers, k->num_kmers, sizeof(char*), idCmp ) !=
	 NULL ) {
      free( test_kmer );
      return 1;
    }
  }

  free( test_kmer );
  return 0;

}

