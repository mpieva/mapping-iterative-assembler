#include "mia.h"



/* reiterate_assembly
   Args: (1) a pointer to a sequence to be used as the new reference
         (2) a MapAlignmentP big enough to store all the alignments
	 (3) a FSDB with sequences to be realigned
	 (4) a AlignmentP big enough for the alignments
	 (5) a front PWAlnFragP for storing front alignments
	 (6) a back PWAlnFragP for storing back alignments
	 (7) a PSSMP with the forward substitution matrices
	 (8) a PSSMP with the revcom substitution matrices
   Aligns all the FragSeqs from fsdb to the new reference, using the
   as and ae fields to narrow down where the alignment happens
   Resets the maln and writes all the results there
   Returns void
*/
void reiterate_assembly( char* new_ref_seq, int iter_num,
			 MapAlignmentP maln,
			 FSDB fsdb, AlignmentP a, 
			 PWAlnFragP front_pwaln,
			 PWAlnFragP back_pwaln, 
			 PSSMP ancsubmat,
			 PSSMP rcancsubmat ) {
  int i, j,
    ref_len,
    ref_start, 
    ref_end,
    ref_frag_len, 
    max_score,
    rc_score,
    aln_seq_len;
  FragSeqP fs;
  char iter_ref_id[MAX_ID_LEN + 1];
  char tmp_rc[INIT_ALN_SEQ_LEN + 1];
  char iter_ref_desc[] = "iteration assembly";

  /* Set up maln->ref
     Keep his seq separate from the external assembly because that
     is malloced and freed elsewhere
  */
  sprintf( iter_ref_id, "ConsAssem.%d", iter_num );
  free( maln->ref->seq );
  if ( maln->ref->rcseq != NULL ) {
    free( maln->ref->rcseq );
  }
  free( maln->ref->gaps );

  ref_len = strlen( new_ref_seq );
  maln->ref->seq = (char*)save_malloc((ref_len + 1)* sizeof(char));
  strcpy( maln->ref->seq, new_ref_seq );
  maln->ref->rcseq = NULL; // never again!
  /* Keep the ID and description the same if this is the 1st
     iteration. Otherwise, set it to the generic ones */
  if ( iter_num > 1 ) {
    strcpy( maln->ref->id, iter_ref_id );
    strcpy( maln->ref->desc, iter_ref_desc );
  }

  maln->ref->seq_len = ref_len;
  maln->ref->size = (ref_len+1);

  if ( maln->ref->circular ) {
    add_ref_wrap( maln->ref );
  }
  else {
    maln->ref->wrap_seq_len = maln->ref->seq_len;
  }
  maln->ref->gaps = 
    (int*)save_malloc((maln->ref->wrap_seq_len+1) * sizeof(int));
  for( i = 0; i <= maln->ref->wrap_seq_len; i++ ) {
    maln->ref->gaps[i] = 0;
  }

  /* Reset its AlnSeqArray ->ins to all point to null */
  for ( i = 0; i < maln->num_aln_seqs; i++ ) {
    aln_seq_len = strlen(maln->AlnSeqArray[i]->seq);
    for ( j = 0; j < aln_seq_len; j++ ) {
      /* We couldn't have malloced any sequence for
	 inserts past our length; anything non-NULL
	 out there is cruft */
      if ( maln->AlnSeqArray[i]->ins[j] != NULL ) {
	free( maln->AlnSeqArray[i]->ins[j] );
	maln->AlnSeqArray[i]->ins[j] = NULL;
      }
    }
  }

  /* Now, remake the hpcl and hprl arrays if hp_special */
  if ( a->hp ) {
    free( a->hpcl );
    free( a->hpcs );
    a->hpcl = (int*)save_malloc(maln->ref->wrap_seq_len*sizeof(int));
    a->hpcs = (int*)save_malloc(maln->ref->wrap_seq_len*sizeof(int));
    pop_hpl_and_hps( maln->ref->seq, 
		     maln->ref->wrap_seq_len,
		     a->hpcl, a->hpcs );     
  }

  /* Reset the number of aligned sequences in the maln */
  maln->num_aln_seqs = 0;

  /* OK, ref is set up. Let's go through all the sequences in fsdb
     and re-align them to the new reference. 
     If it's a revcom alignment,
     just use the rcancsubmat */
  for( i = 0; i < fsdb->num_fss; i++ ) {
    fs = fsdb->fss[i];

    /* Special case of distant reference and 
       !fs->strand_known => try to realign both strands
       against the entire reference to learn the 
       strand and alignment region
    */
    if ( maln->distant_ref &&
	 (fs->strand_known == 0 ) &&
	 (iter_num > 1) ) {
      ref_start = 0;
      ref_end = maln->ref->wrap_seq_len;
      ref_frag_len = ref_end - ref_start;
      a->seq1 = &maln->ref->seq[0];
      a->len1 = ref_frag_len;
      pop_s1c_in_a( a );
      a->seq2 = fs->seq;
      a->len2 = strlen( a->seq2 );
      pop_s2c_in_a( a );
      if ( a->hp ) {
	pop_hpl_and_hps( a->seq2, a->len2, a->hprl, a->hprs );
	pop_hpl_and_hps( a->seq1, a->len1, a->hpcl, a->hpcs );
      }
      /* Align it! */
      dyn_prog( a );
      /* Find the best forward score */
      max_score = max_sg_score( a );
      if ( max_score > FIRST_ROUND_SCORE_CUTOFF ) {
	fs->strand_known = 1;
	fs->rc = 0;
	find_align_begin( a );
	fs->as = a->abc;
	fs->ae = a->aec;
	fs->score = max_score;
      }

      /* Now, try reverse complement */
      aln_seq_len = strlen( fs->seq );
      a->submat = rcancsubmat;
      for ( j = 0; j < aln_seq_len; j++ ) {
	tmp_rc[j] = revcom_char(fs->seq[aln_seq_len-(j+1)]);
      }
      tmp_rc[aln_seq_len] = '\0';
      a->seq2 = tmp_rc;
      pop_s2c_in_a( a );
      if ( a->hp ) {
	pop_hpl_and_hps( a->seq2, a->len2, a->hprl, a->hprs );
	pop_hpl_and_hps( a->seq1, a->len1, a->hpcl, a->hpcs );
      }
      dyn_prog( a );
      max_score = max_sg_score( a );
      if ( (max_score > FIRST_ROUND_SCORE_CUTOFF) &&
	   (max_score > fs->score) ) {
	fs->strand_known = 1;
	fs->rc = 1;
	find_align_begin( a );
	fs->as = a->abc;
	fs->ae = a->aec;
	fs->score = max_score;
	strcpy( fs->seq, tmp_rc );
      }
    }

    /* Do we know the strand (either because we've always
       known it or we just learned it, doesn't matter) */
    if ( fs->strand_known ) {
      if ( fs->rc ) {
	a->submat = rcancsubmat;
      }
      else {
	a->submat = ancsubmat;
      }

      a->seq2 = fs->seq;
      a->len2 = strlen( a->seq2 );
      pop_s2c_in_a( a );

      /* Set up the alignment limits on the reference */
      if ( ((fs->as - REALIGN_BUFFER) < 0 ) ) {
	ref_start = 0;
      }
      else {
	ref_start = (fs->as - REALIGN_BUFFER);
      }
      if ( (fs->ae + REALIGN_BUFFER + 1) > 
	   maln->ref->wrap_seq_len ) {
	ref_end = maln->ref->wrap_seq_len;
      }
      else {
	ref_end = fs->ae + REALIGN_BUFFER;
      }

      /* Check to make sure the regions encompassed by ref_start to
	 ref_end is reasonable given how long this fragment is. If
	 not, just realign this whole mofo again because the reference
	 has probably changed a lot between iterations */
      if ( (ref_start + a->len2) > ref_end ) {
	ref_start = 0;
	ref_end = maln->ref->wrap_seq_len;
      }
    
      ref_frag_len = ref_end - ref_start;
      a->seq1 = &maln->ref->seq[ref_start];
      a->len1 = ref_frag_len;
      pop_s1c_in_a( a );
      
      /* If we want the homopolymer discount, the necessary arrays of
	 hp starts and lengths must be set up anew */
      if ( a->hp ) {
	pop_hpl_and_hps( a->seq2, a->len2, a->hprl, a->hprs );
	pop_hpl_and_hps( a->seq1, a->len1, a->hpcl, a->hpcs );
      }

      /* Align it! */
      dyn_prog( a );
    
      /* Find the best score */
      max_score = max_sg_score( a );

      find_align_begin( a );

      /* First, put all alignment in front_pwaln */
      populate_pwaln_to_begin( a, front_pwaln );
      
      /* Load up front_pwaln */
      strcpy( front_pwaln->ref_id, maln->ref->id );
      strcpy( front_pwaln->ref_desc, maln->ref->desc );
      
      strcpy( front_pwaln->frag_id, fs->id );
      strcpy( front_pwaln->frag_desc, fs->desc );
      
      front_pwaln->trimmed = fs->trimmed;
      front_pwaln->revcom  = fs->rc;
      front_pwaln->num_inputs = fs->num_inputs;
      front_pwaln->segment = 'a';
      front_pwaln->score = a->best_score;
  
      front_pwaln->start = a->abc + ref_start;
      front_pwaln->end   = a->aec + ref_start;

      /* Update stats for this FragSeq */
      fs->as = a->abc + ref_start;
      fs->ae = a->aec + ref_start;
      fs->unique_best = 1;
      fs->score = a->best_score;

      if ( front_pwaln->end > maln->ref->seq_len ) {
	/* This alignment wraps around - adjust the end to
	   demonstrate this for split_maln check */
	front_pwaln->end = front_pwaln->end - maln->ref->seq_len;
      }

      if ( front_pwaln->start > front_pwaln->end ) {
	/* Move wrapped bit to back_pwaln */
	split_pwaln( front_pwaln, back_pwaln, maln->ref->seq_len );
	merge_pwaln_into_maln( front_pwaln, maln );
	fs->front_asp = maln->AlnSeqArray[maln->num_aln_seqs - 1];
	merge_pwaln_into_maln( back_pwaln, maln );
	fs->back_asp = maln->AlnSeqArray[maln->num_aln_seqs - 1];
      }
      else { 
	merge_pwaln_into_maln( front_pwaln, maln );
	fs->front_asp = maln->AlnSeqArray[maln->num_aln_seqs - 1];
      }
    }
  }
  return;
}


/* all_lower
   Args: (1) Pointer to char array (seq)
         (2) int number of characters to check (len)
   Returns: int 1 => first len charaters in seq are all lower case
                0 => at least one of the characters is not lower case
*/
inline int all_lower( const char* seq, const int kmer_len ) {
  size_t i;
  for( i = 0; i < kmer_len; i++ ) {
    if ( isupper( seq[i] ) ) {
      return 0;
    }
  }
  return 1;
}


void help( void ) {
  printf( "\n\n%s -- Mapping Iterativ Assembler V %s\n",PACKAGE_NAME, PACKAGE_VERSION);
  printf( "       A tool for creating short read assemblies.\n\n");
  printf( "Copyright Richard E. Green, Michael Siebauer 2008-2009\n");
  printf( "Report bugs to <%s>.\n",PACKAGE_BUGREPORT);
  printf( "===============================+++++++++++++==\n");
  printf( "\nUsage:\n");
  printf( "mia -r <reference sequence>\n" );
  printf( "    -f <fasta or fastq file of fragments to align>\n" );
  printf( "    -s <substitution matrix file> (if not supplied an default matrix is used)\n" );
  printf( "    -m <root file name for maln output file(s)> (assembly.maln.iter)\n" );
  printf( "    \nFILTER parameters:\n" );
  printf( "    -u fasta database has repeat sequences, keep one based on alignment score\n" );
  printf( "    -U fasta database has repeat sequences, keep one based on sum of q-scores\n" );
  printf( "    -C collapse sequences with same start, end, strand info into a single sequence\n" );
  printf( "    -A use adapter presence and coordinate information to more aggressively\n" );
  printf( "       remove repeat sequences - suitable only for 454 sequences that have not\n" );
  printf( "       already been adapter trimmed\n" );
  printf( "    -T fasta database has adapters, trim these\n" );
  printf( "    -a <adapter sequence or code>\n" );
  printf( "    -k <use kmer filter with kmers of this length>\n" );
  printf( "    -I <filename of list of sequence IDs to use, ignoring all others>\n" );
  printf( "    \nALIGNMENT parameters:\n" );
  printf( "    -p <consensus calling code; default = 1>\n" );
  printf( "    -c means reference/assembly is circular\n" );
  printf( "    -i iterate assembly until convergence\n" );
  printf( "    -F <only output the FINAL assembly, not each iteration>\n" );
  printf( "    -D <distantly related reference sequence>\n" );
  printf( "    -h give special discount for homopolymer gaps\n" );
  printf( "    -M <use lower-case soft-masking of kmers>\n" );
  printf( "    -H <do not do dynamic score cutoff, instead use this Hard score cutoff>\n" );
  printf( "    -S <slope of length/score cutoff line>\n" );
  printf( "    -N <intercept of length/score cutoff line>\n" );
  printf( "The default substitution matrix used the following parameters:\n" );
  printf( "  MATCH=%d, MISMATCH=%d, N=%d for all positions\n", FLAT_MATCH, FLAT_MISMATCH, N_SCORE);

  printf( "The procedure for removing bad-scoring alignments from the assembly is:\n" );
  printf( "Default: fit a line to length versus score and remove reads that are\n" );
  printf( "less that SCORE_CUTOFF_BUFFER than the average score for its length.\n" );
  printf( "If -H is specified then this hard score cutoff is applied to all reads.\n" );
  printf( "This is preferable if all reads are the same length.\n" );
  printf( "If -S or -N are specified, then these are used as the slope and intercept\n" );
  printf( "of a length/score line. Reads must score above this line to be included.\n" );
  printf( "If only one of -S or -N is specified then the default values are used for\n" );
  printf( "the other (default S = %.1f; default N = %.1f)\n", DEF_S, DEF_N );
  printf( "The kmer filter requires that a sequence fragment have at least one\n" );
  printf( "kmer of the specified length in common with the reference sequence in\n" );
  printf( "order to align it. For 36nt Solexa data, a value of 12 works well.\n" );
  printf( "The -p option specifies how the new consensus assembly sequence is called\n" );
  printf( "at each iteration:\n" );
  printf( "1 => Any base whose aggregate score is MIN_SC_DIFF_CONS better than all\n" );
  printf("      others is the assembly base. If none is, then N is the assembly base.\n" );
  printf( "2 => The best scoring base whose aggregate score is better than MIN_SCORE_CONS\n" );
  printf( "     is the assembly base. If none is, then N is the assembly base.\n" );
  printf( "If -T is specified, mia will attempt to find and trim adapters on\n" );
  printf( "each sequence. The adapter sequence itself can be specified by a\n" );
  printf( "one letter code as argument to -a. N or n => Neandertal adapter\n" );
  printf( "                  any other single letter => Standard GS FLX adapter\n" );
  printf( "              sequence (less than 127 nt) => user-specified adapter\n" );
}

int main( int argc, char* argv[] ) {

  char mat_fn[MAX_FN_LEN+1];
  char maln_fn[MAX_FN_LEN+1];
  char fastq_out_fn[MAX_FN_LEN+1];
  char maln_root[MAX_FN_LEN+1];
  char ref_fn[MAX_FN_LEN+1];
  char frag_fn[MAX_FN_LEN+1];
  char adapter_code[2]; // place to keep the argument for -a (which adapter to trim)
  char* c_time; // place to keep asctime string
  char* test_id;

  int ich;
  int any_arg = 0;
  int Hard_cut = 0; // If 0 => use dynamic score cutoff, if > 0 use this instead
  int circular = 0; // Boolean, TRUE if reference sequence is circular
  int make_fastq = 0; // Boolean, TRUE if we should also output fastq database of seqs in assembly
  int seq_code = 0; // code to indicate sequence input format; 0 => fasta; 1 => fastq
  int do_adapter_trimming = 0; // Boolean, TRUE if we should try to trim
                               // adapter from input sequences
  int iterate = 0; //Boolean, TRUE means interate the assembly until convergence
  // on an assembled sequence
  int FINAL_ONLY = 0; //Boolean, TRUE means only write out the final assembly maln file
                      //         FALSE (default) means write out each one
  int ids_rest = 0; // Boolean, TRUE means restrict analysis to IDs in input file
  int repeat_filt = 0; //Boolean, TRUE means remove sequences that are repeats, 
                       // keeping best-scoring representative
  int repeat_qual_filt = 0; //Boolean, TRUE means remove sequences that are repeats,
                            // keeping best quality score sum representative
  int just_outer_coords = 1; // Boolean, TRUE means just use strand, start, and end to
                             // determine if sequences are redundant
  int SCORE_CUT_SET = 0; //Boolean, TRUE means user has set a length/score cutoff line
  int seen_seqs = 0;
  int hp_special = 0; // Boolean, TRUE means user wants hp gap special discount
  int distant_ref = 0; // Boolean, TRUE means the initial reference sequence is
                       // known to be distantly related so keep trying to align all
                       // sequences each round
  int kmer_filt_len = -1; // length of kmer filtering, if user wants it; otherwise
                          // special value of -1 indicates this is unset
  int soft_mask = 0; //Boolean; TRUE => do not use kmers that are all lower-case
                     //        FALSE => DO use all kmers, regardless of case
  int iter_num; // Number of iterations of assembly done
  int collapse = 0; // Boolean; TRUE => collapse input sequences in FSDB to improve
                    //                  sequence quality
                    //          FALSE => (default) keep all sequences
  double slope     = DEF_S; // Set these to default unless, until user changes
  double intercept = DEF_N; // them 
  MapAlignmentP maln, // Contains all fragments initially better
                      // than FIRST_ROUND_SCORE_CUTOFF
    culled_maln;      // Contains all fragments with scores
                      // better than SCORE_CUTOFF
  AlignmentP fw_align, rc_align, adapt_align;
  
  PSSMP ancsubmat   = init_flatsubmat();
  PSSMP rcancsubmat = revcom_submat(ancsubmat);
  const PSSMP flatsubmat  = init_flatsubmat();

  KPL* fkpa; // Place to keep forward kmer array if user requested kmer 
  KPL* rkpa; // Place to keep reverse kmer array if user requested kmer 
  IDsListP good_ids;
  FragSeqP frag_seq;
  PWAlnFragP front_pwaln, back_pwaln;
  FSDB fsdb; // Database to hold sequences to iterate over
  FILE* FF;
  time_t curr_time;


  char maln_root_def[] = "assembly.maln.iter";
  extern int optind;
  extern char* optarg;
  char neand_adapt[] = "GTCAGACACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGG";
  char stand_adapt[] = "CTGAGACACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGG";
  char user_def_adapt[128];
  char* adapter; // set to either neand_adapt or stand_adapt based on user preference
  adapter = neand_adapt; // Default is Neandertal
  char* assembly_cons;
  char* last_assembly_cons;
  int cc = 1; // consensus code for calling consensus base
  int i;

  /* Set the default output filename until the user overrides it */
  strcpy( maln_root, maln_root_def );


  /* Process command line arguments */
  while( (ich=getopt( argc, argv, "s:r:f:m:a:p:H:I:S:N:k:q:FTciuhDMUAC" )) != -1 ) {
    switch(ich) {
    case 'c' :
      circular = 1;
      break;
    case 'q' :
      make_fastq = 1;
      strcpy( fastq_out_fn, optarg );
    case 'C' :
      collapse = 1;
      break;
    case 'i' :
      iterate = 1;
      break;
    case 'h' :
      hp_special = 1;
      break;
    case 'u' :
      repeat_filt = 1;
      break;
    case 'A' :
      just_outer_coords = 0;
      break;
    case 'U' :
      repeat_qual_filt = 1;
      break;
    case 'D' :
      distant_ref = 1;
      break;
    case 'p' :
      cc = atoi( optarg );
      any_arg = 1;
      break;
    case 'I' :
      good_ids = parse_ids( optarg );
      ids_rest = 1;
      break;
    case 'H' :
      Hard_cut = atoi( optarg );
      if ( Hard_cut <= 0 ) {
	fprintf( stderr, "Hard cutoff (-H) must be positive\n" );
	help();
	exit( 0 );
      }
      any_arg = 1;
      break;
    case 'M' :
      soft_mask = 1;
      break;
    case 's' :
      strcpy( mat_fn, optarg );
      free( ancsubmat ); // trash the flat submat we initialized with
      ancsubmat   = read_pssm( mat_fn );
      free( rcancsubmat ); // trash the init rcsubmat, too
      rcancsubmat = revcom_submat( ancsubmat );
      any_arg = 1;
      break;
    case 'r' :
      strcpy( ref_fn, optarg );
      any_arg = 1;
      break;
    case 'k' :
      kmer_filt_len = atoi( optarg );
      any_arg = 1;
      break;
    case 'f' :
      strcpy( frag_fn, optarg );
      any_arg = 1;
      break;
    case 'm' :
      strcpy( maln_root, optarg );
      any_arg = 1;
      break;
    case 'T' :
      do_adapter_trimming = 1;
      break;
    case 'a' :
      if ( strlen( optarg ) > 127 ) {
	  fprintf( stderr, "That adapter is too big!\nMIA will use the standard adapter.\n" );
	  adapter = stand_adapt;
      }
      else {
	strcpy( user_def_adapt, optarg );
	  if ( strlen( user_def_adapt ) > 1 ) {
	    adapter = user_def_adapt;
	  }
	  else {
	    if ( !( (user_def_adapt[0] == 'n') ||
		    (user_def_adapt[0] == 'N') ) ) {
	      adapter = stand_adapt;
	    }
	    else {
	      adapter = neand_adapt;
	    }
	  }
      }
      break;
    case 'S' :
      slope = atof( optarg );
      SCORE_CUT_SET = 1;
      break;
    case 'N' :
      intercept = atof( optarg );
      SCORE_CUT_SET = 1;
      break;
    case 'F' :
      FINAL_ONLY = 1;
      break;
    default :
      help();
      exit( 0 );
    }
  }

  if ( !any_arg ) {
    help();
    exit( 0 );
  }

  if ( optind != argc ) {
    fprintf( stderr, "There seems to be some extra cruff on the command line that mia does not understand.\n" );
  }

  /* Start the clock... */
  curr_time = time(NULL);
  //  c_time = (char*)save_malloc(64*sizeof(char));
  //  c_time = asctime(localtime(&curr_time));

  /* Announce that we're starting */
  fprintf( stderr, 
	   "Starting assembly of %s\nusing %s\nas reference at %s\n", 
	   frag_fn, ref_fn, 
	   asctime(localtime(&curr_time)) );


  /* Set up the maln structure */
  maln = (MapAlignmentP)init_map_alignment();
  maln->cons_code = cc; 
  if ( maln == NULL ) {
    fprintf( stderr, "Not enough memories for this\n" );
    exit( 1 );
  }

  /* Set the distant_ref flag */
  maln->distant_ref = distant_ref;

  /* Set up the FSDB for keeping good-scoring sequence in memory */
  fsdb = init_FSDB();
  if ( fsdb == NULL ) {
    fprintf( stderr, "Not enough memories for holding sequences\n" );
    exit( 1 );
  }

  /* Read in the reference sequence and make reverse complement, too*/
  if ( read_fasta_ref( maln->ref, ref_fn ) != 1 ) {
    fprintf( stderr, "Problem reading reference sequence file %s\n", ref_fn );
    exit( 1 );
  }

  /* Add wrap-around sequence (rc, too) and set maln->ref->circular
     if it's circular */
  if ( circular ) {
    add_ref_wrap( maln->ref );
  }
  else {
    maln->ref->wrap_seq_len = maln->ref->seq_len;
  }
  /* Add space for the gaps array */
  maln->ref->gaps = (int*)save_malloc((maln->ref->wrap_seq_len+1) *
				      sizeof(int));
  for( i = 0; i <= maln->ref->wrap_seq_len; i++ ) {
    maln->ref->gaps[i] = 0;
  }

  /* Set up fkpa and rkpa for list of kmers in the reference (forward and
     revcom strand) if user wants kmer filtering */
  if ( kmer_filt_len > 0 ) {
    fprintf( stderr, "Making kmer list for k-mer filtering...\n" );
    fkpa = init_kpa(kmer_filt_len);
    rkpa = init_kpa(kmer_filt_len);
    /* 
    kmer_list = (KmersP)pop_kmers( maln->ref, kmer_filt_len );
    */
    populate_kpa( fkpa, maln->ref->seq, 
		  maln->ref->wrap_seq_len, kmer_filt_len, 
		  soft_mask );
    populate_kpa( rkpa, maln->ref->rcseq, 
		  maln->ref->wrap_seq_len, kmer_filt_len,
		  soft_mask );
  }

  /* Now kmer arrays have been made if requested. We can upper case
     the reference sequences. */
  make_ref_upper( maln->ref );

  /* Set up FragSeqP to point to a FragSeq */
  frag_seq = (FragSeqP)save_malloc(sizeof(FragSeq));

  /* Set up the alignment structures for forward and reverse
     complement alignments */
  fw_align = (AlignmentP)init_alignment( INIT_ALN_SEQ_LEN,
					 (maln->ref->wrap_seq_len + 
					  (2*INIT_ALN_SEQ_LEN)),
					 0, hp_special );
  rc_align = (AlignmentP)init_alignment( INIT_ALN_SEQ_LEN,
					 (maln->ref->wrap_seq_len + 
					  (2*INIT_ALN_SEQ_LEN)),
					 1, hp_special );

  /* Set up the alignment structure for adapter trimming, if user
     wants that */
  if ( do_adapter_trimming ) {
    adapt_align = (AlignmentP)init_alignment( INIT_ALN_SEQ_LEN,
					      INIT_ALN_SEQ_LEN,
					      0, hp_special );
    /* Setup the flatsubmat */
    //flatsubmat = init_flatsubmat();
    adapt_align->submat = flatsubmat;

    adapt_align->seq2   = adapter;
    adapt_align->len2   = strlen( adapt_align->seq2 );
    pop_s2c_in_a( adapt_align );
    if ( hp_special ) {
      pop_hpl_and_hps( adapt_align->seq2, adapt_align->len2,
		       adapt_align->hprl, adapt_align->hprs );
    }
    /* Set for a semi-global that pays a penalty for unaligning the
       beginning of the adapter, but not for the end of the adapter.
       This is because if the sequence read (align->seq1) ends, then
       we won't see any more of the adapter. When we search for the
       best alignment, we'll only look in the last column, requiring that
       all of align->seq1 is accounted for */
    adapt_align->sg5    = 1;
    adapt_align->sg3    = 0;
  }

  fw_align->seq1 = maln->ref->seq;
  rc_align->seq1 = maln->ref->rcseq;
  if ( circular ) {
    fw_align->len1 = maln->ref->wrap_seq_len;
    rc_align->len1 = maln->ref->wrap_seq_len;
  }
  else {
    fw_align->len1 = maln->ref->seq_len;
    rc_align->len1 = maln->ref->seq_len;
  }

  /* Now the reference sequence and its reverse complement are
     prepared, put the s1c lookup codes in */
  pop_s1c_in_a( fw_align );
  pop_s1c_in_a( rc_align );

  if ( hp_special ) {
    pop_hpl_and_hps( fw_align->seq1, fw_align->len1,
		     fw_align->hpcl, fw_align->hpcs );
    pop_hpl_and_hps( rc_align->seq1, rc_align->len1,
		     rc_align->hpcl, rc_align->hpcs );
  }

  /* One by one, go through the input file of fragments to be aligned.
     Align them to the reference. For each fragment generating an
     alignment score better than the cutoff, merge it into the maln
     alignment. Keep track of those that don't, too. */
  FF = fileOpen( frag_fn, "r" );
  seq_code = find_input_type( FF );

  //LOG = fileOpen( log_fn, "w" );
  front_pwaln = (PWAlnFragP)save_malloc( sizeof(PWAlnFrag));
  back_pwaln  = (PWAlnFragP)save_malloc( sizeof(PWAlnFrag));

  /* Give some space to remember the IDs as we see them */
  test_id = (char*)save_malloc(MAX_ID_LEN * sizeof(char));

  /* Announce we're strarting alignment of fragments */
  fprintf( stderr, "Starting to align sequences to the reference...\n" );

  while( read_next_seq( FF, frag_seq, seq_code ) ) {
    seen_seqs++;
    strcpy( test_id, frag_seq->id );
    if ( DEBUG ) {
      fprintf( stderr, "%s\n", frag_seq->id );
    }
    if ( !ids_rest ||
	 ( bsearch( &test_id, good_ids->ids, 
		    good_ids->num_ids,
		    sizeof(char*), idCmp ) 
	   != NULL ) ) {

      if ( do_adapter_trimming ) {
	/* Trim sequence (set frag_seg->trimmed and 
	   frag_seg->trim_point field) */
	trim_frag( frag_seq, adapter, adapt_align );
      }
      else {
	frag_seq->trimmed = 0;
      }

      /* Check if kmer filtering. If so, filter */
      if ( new_kmer_filter( frag_seq, fkpa, rkpa, kmer_filt_len,
			    fw_align, rc_align ) ) {
	/* Align this fragment to the reference and write 
	   the result into pwaln; use the ancsubmat, not the reverse
	   complemented rcsancsubmat during this first iteration because
	   all sequence is forward strand
	*/
	fw_align->submat = ancsubmat;
	rc_align->submat = ancsubmat;
	
	if ( sg_align( maln, frag_seq, fsdb, 
		       fw_align, rc_align,
		       front_pwaln, 
		       back_pwaln ) == 0 ) {
	  fprintf( stderr, "Problem handling %s\n", frag_seq->id );
	}
      }  
    }
    if ( seen_seqs % 1000 == 0 ) {
      fprintf( stderr, "." );
    }
    if ( seen_seqs % 80000 == 0 ) {
      fprintf( stderr, "\n" );
    }
  }

  /* Now, fsdb is complete and points to all the things in maln.
     So we can fill in the AlnSeqP->smp array for everything in the 
     maln->AlnSeqArray to know which matrices to use for *CALLING* 
     a consensus; Conveniently, there are pointers to all of these
     in the fss->fss[X]->front|back_asp */
  pop_smp_from_FSDB( fsdb, PSSM_DEPTH );

  //fprintf( LOG, "__Finished with initial alignments__" );
  //fflush( LOG );
  fprintf( stderr, "\n" );
  iter_num = 1;

  /* Now, we need a new MapAlignment, culled_maln, that is big
     enough to hold all the unique guys from maln */
  culled_maln = init_culled_map_alignment( maln );

  /* Filtering repeats announcement */
  fprintf( stderr, "Repeat and score filtering\n" );

  /* If user wants to filter against repeats by alignment score, do it */
  if ( repeat_filt ) {  
    /* Sort fsdb by fsdb->as */
    sort_fsdb( fsdb );
    
    /* Now, everything is sorted in fsdb, so I can easily see
       which guys are unique by as, ae, and rc fields */
    set_uniq_in_fsdb( fsdb, just_outer_coords );
  }

  /* If user wants to filter against repeats by q-score sum, do it */
  if ( repeat_qual_filt ) {  
    /* Sort fsdb by fsdb->as */
    sort_fsdb_qscore( fsdb );
    
    /* Now, everything is sorted in fsdb, so I can easily see
       which guys are unique by as, ae, and rc fields */
    set_uniq_in_fsdb( fsdb, just_outer_coords );
  }

  /* Now, we know which sequences are unique, so make a
     culled_maln with just the unique guys */
  cull_maln_from_fsdb( culled_maln, fsdb, Hard_cut, 
		       SCORE_CUT_SET, slope, intercept );

  fclose(FF);

  /* Tell the culled_maln which matrices to use for assembly */
  culled_maln->fpsm = ancsubmat;
  culled_maln->rpsm = rcancsubmat;

  sort_aln_frags( culled_maln ); //invalidates fsdb->front|back_asp fields!

  fw_align->submat = ancsubmat;
  fw_align->sg5 = 1;
  fw_align->sg3 = 1;

  last_assembly_cons = (char*)save_malloc((maln->ref->seq_len +1) * 
				     sizeof(char));
  strncpy( last_assembly_cons, maln->ref->seq, 
	   maln->ref->seq_len );
  last_assembly_cons[maln->ref->seq_len] = '\0';

  /* Re-align everything with revcomped
     sequence and substitution matrices, but first
     unmask all alignment positions and collapse sequences
     if requested
  */
  memset(fw_align->align_mask, 1, fw_align->len1);
  if ( collapse ) {
    collapse_FSDB( fsdb, Hard_cut, SCORE_CUT_SET, 
		   slope, intercept );
  }
  reiterate_assembly( last_assembly_cons, iter_num, maln, fsdb,
		      fw_align, front_pwaln, back_pwaln, 
		      ancsubmat, rcancsubmat );
  pop_smp_from_FSDB( fsdb, PSSM_DEPTH );
  fprintf( stderr, "Repeat and score filtering\n" );
  if ( repeat_filt ) {
    sort_fsdb( fsdb );
    set_uniq_in_fsdb( fsdb, just_outer_coords );
  }
  if ( repeat_qual_filt ) {  
    sort_fsdb_qscore( fsdb );
    set_uniq_in_fsdb( fsdb, just_outer_coords );
  }
  cull_maln_from_fsdb( culled_maln, fsdb, Hard_cut,
		       SCORE_CUT_SET, slope, intercept );
  
  
  /* Tell the culled_maln which matrices to use for assembly */
  culled_maln->fpsm = ancsubmat;
  culled_maln->rpsm = rcancsubmat;
  
  //invalidates fsdb->front|back_asp fields!
  sort_aln_frags( culled_maln );
  sprintf( maln_fn, "%s.%d", maln_root, iter_num );
  if ( !iterate || !FINAL_ONLY ) {
    write_ma( maln_fn, culled_maln );
    if ( make_fastq ) {
      write_fastq( fastq_out_fn, fsdb );
    }
  }

  /* Are we iterating (re-aligning to the a new consensus? */
  if (iterate) {
    /* New assembly consensus announcement */
    fprintf( stderr, "Generating new assembly consensus\n" );
    assembly_cons = consensus_assembly_string( culled_maln );

    while( ( strcmp( assembly_cons, last_assembly_cons ) != 0) &&
	   (iter_num < MAX_ITER) ) {
      /* Another round...*/
      iter_num++;
      free( last_assembly_cons );
      last_assembly_cons = assembly_cons;

      fprintf( stderr, "Starting assembly iteration %d\n", 
	       iter_num );

      /* If the user wants collapsed sequences, now is the time */
      if ( collapse ) {
	collapse_FSDB( fsdb, Hard_cut, SCORE_CUT_SET, 
		       slope, intercept );
      }

      reiterate_assembly( assembly_cons, iter_num, maln, fsdb, 
			  fw_align, front_pwaln, back_pwaln,
			  ancsubmat, rcancsubmat );

      pop_smp_from_FSDB( fsdb, PSSM_DEPTH );

      fprintf( stderr, "Repeat and score filtering\n" );
      if ( repeat_filt ) {
	sort_fsdb( fsdb );
	set_uniq_in_fsdb( fsdb, just_outer_coords );
      }
      if ( repeat_qual_filt ) {
	sort_fsdb_qscore( fsdb );
	set_uniq_in_fsdb( fsdb, just_outer_coords );
      }
      cull_maln_from_fsdb( culled_maln, fsdb, Hard_cut,
			   SCORE_CUT_SET, slope, intercept );

      
      /* Tell the culled_maln which matrices to use for assembly */
      culled_maln->fpsm = ancsubmat;
      culled_maln->rpsm = rcancsubmat;

      //invalidates fsdb->front|back_asp fields!
      sort_aln_frags( culled_maln );

      sprintf( maln_fn, "%s.%d", maln_root, iter_num );
      if ( !FINAL_ONLY ) {
	fprintf( stderr, "Writing maln file for iteration %d\n", 
		 iter_num );
	write_ma( maln_fn, culled_maln );
      }
      assembly_cons = consensus_assembly_string( culled_maln );
    }
  
    /* Convergence? */
    if ( strcmp( assembly_cons, last_assembly_cons ) == 0 ) {
      fprintf( stderr, "Assembly convergence - writing final maln\n" );
      write_ma( maln_fn, culled_maln );
    }
    else {
      fprintf( stderr, "Assembly did not converge after % rounds, quitting\n" );
      write_ma( maln_fn, culled_maln );
    }
    if ( make_fastq ) {
      write_fastq( fastq_out_fn, fsdb );
    }
  }

  /* No iteration, but we must still re-align everything with revcomped
     sequence and substitution matrices to keep scores comparable to what
     they would have been had we iterated */

  /* Announce we're finished */
  curr_time = time(NULL);
  //  c_time    = asctime(localtime(&curr_time));
  fprintf( stderr, "Assembly finished at %s\n",
	   asctime(localtime(&curr_time)) );

  exit( 0 );
}
