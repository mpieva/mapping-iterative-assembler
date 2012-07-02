#include "io.h"

/* find_input_type
   Args: 1. FILE* pointer to file to be analyzed
   Returns: sequence code indicating what kind of sequence file
            this is:
	    0 => fasta
	    1 => fastq
   Resets the input FILE pointer to the beginning of the file
*/
int find_input_type( FILE * FF ) {
  char c;
  c = fgetc( FF );
  ungetc( c, FF );  
  if ( c == '@' ) {
    return 1;
  }

  if ( c == '>' ) {
    return 0;
  }

  /* default */
  return 0;
}


/* read_next_seq
   Args: 1. FILE* pointer to file being read
         2. FragSeqP pointer to FragSeq where the next sequence data will go
	 3. int code indicating which parser to use
   Returns: TRUE if a sequence was read,
            FALSE if EOF
*/
int read_next_seq( FILE * FF, FragSeqP frag_seq, int seq_code ) {
  if ( seq_code == 0 ) {
    return read_fasta( FF, frag_seq );
  }
  if ( seq_code == 1 ) {
    return read_fastq( FF, frag_seq );
  }
}

/* read_fastq
   Args 1. pointer to file to be read
        2. pointer to FragSeq to put the sequence into
   Returns: TRUE if a sequence was read,
            FALSE if EOF
*/
int read_fastq ( FILE * fastq, FragSeqP frag_seq ) {
  char c;
  size_t i;
  c = fgetc( fastq );
  if ( c == EOF ) return 0;
  if ( c != '@' ) {
    fprintf( stderr, "While reading fastq file, saw record not beginning with @\n" );
    fprintf( stderr, "Maybe badly formed input? Continuing, anyway...\n" );
    return 0;
  }

  /* get identifier */
  i = 0;
  while( (!isspace(c=fgetc( fastq ) ) &&
	  (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return 0;
    }
    frag_seq->id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      /* Id is too long - truncate it now */
      frag_seq->id[i] = '\0';
    }
  }
  frag_seq->id[i] = '\0';

  /* Now, everything else on the line is description (if anything)
     although fastq does not appear to formally support description */
  if ( c == '\n' ) {
    frag_seq->desc[0] = '\0';
  }
  else { // some description, uh oh
    while ( (c != '\n') &&
	    (isspace(c)) ) {
      c = fgetc( fastq );
    }
    i = 0;
    while( (c != '\n') &&
	   (i < MAX_DESC_LEN) ) {
      frag_seq->desc[i++] = c;
      c = fgetc( fastq );
    }
    frag_seq->desc[i] = '\0';
  }

  /* Now, read the sequence. This should all be on a single line */
  i = 0;
  c = fgetc( fastq );
  while ( (c != '\n') &&
	  (c != EOF) &&
	  (i < INIT_ALN_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      frag_seq->seq[i++] = c;
    }
    c = fgetc( fastq );
  }
  frag_seq->seq[i] = '\0';
  frag_seq->seq_len = i;
  /* If the reading stopped because the sequence was longer than
     INIT_ALN_SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == INIT_ALN_SEQ_LEN ) {
    while ( (c != '\n') &&
	    (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  /* Now, read the quality score header */
  c = fgetc( fastq );
  if ( c != '+' ) {
    fprintf( stderr, "Problem reading quality line for %s\n", frag_seq->id );
    return 1;
  }
  /* Zip through the rest of the line, it should be the same identifier
     as before or blank */
  c = fgetc( fastq );
  while( (c != '\n') &&
	 (c != EOF) ) {
    c = fgetc( fastq );
  }

  /* Now, get the quality score line */
  c = fgetc( fastq );
  i = 0;
  while( (c != '\n') &&
	 (c != EOF) &&
	 (i < INIT_ALN_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      frag_seq->qual[i++] = c;
    }
    c = fgetc( fastq );
  }
  frag_seq->qual[i] = '\0';

  frag_seq->qual_sum = calc_qual_sum( frag_seq->qual );

  /* If the reading stopped because the sequence was longer than
     INIT_ALN_SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == INIT_ALN_SEQ_LEN ) {
    while ( (c != '\n') &&
	    (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  if ( i != frag_seq->seq_len ) {
    fprintf( stderr, "%s has unequal sequence and qual line lengths\n", 
	     frag_seq->id );
    return 0;
  }
  return 1;
}

/* calc_qual_sum
   Args: 1. pointer to a string of quality scores for this sequence
   Returns: 1. int - the sum of quality scores for this sequence
   This assumes that quality scores are represented as the 
   ASCII code + 64
*/
inline int calc_qual_sum( const char* qual_str ) {
  size_t i, len;
  int qual_sum = 0;

  len = strlen( qual_str );
  for( i = 0; i < len; i++ ) {
    qual_sum += (qual_str[i] - 33);
  }
  
  return qual_sum;
}


/* read_fasta
   args 1. pointer to file to be read
        2. pointer to FragSeq to put the sequence
   returns: TRUE if sequence was read,
            FALSE if EOF or not fasta
*/
int read_fasta ( FILE * fasta, FragSeqP frag_seq ) {
  char c;
  size_t i;
  c = fgetc( fasta );
  if ( c == EOF ) return 0;
  if ( c != '>' ) return 0;

  /* No quality scores, so initialize this to keep
     stupid valgrind from stupid complaining */
  frag_seq->qual[0] = '\0';

  // get id
  i = 0;
  while( (!isspace( c=fgetc( fasta ) ) &&
	  (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return 0;
    }
    frag_seq->id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      //id is too long - truncate it
      frag_seq->id[i] = '\0';
    }
  }
  frag_seq->id[i] = '\0';

  // everything else on this line is description, if there is anything
  if ( c == '\n' ) {
    frag_seq->desc[0] = '\0';
  }
  else { // not end of line, so skip past the stupid whitespace...
    while( (c != '\n') &&
	   (isspace(c)) ) {
      c = fgetc( fasta );
    }
    ///...everthing else is description
    i = 0;
    ungetc( c, fasta );
    while ( (c != '\n') &&
	    (i < MAX_DESC_LEN) ) {
      frag_seq->desc[i++] = c;
      c = fgetc( fasta );
    }
    frag_seq->desc[i] = '\0';
  }

  // read sequence
  i = 0;
  c = fgetc( fasta );
  while ( ( c != '>' ) &&
          ( c != EOF ) &&
	  (i < INIT_ALN_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      frag_seq->seq[i++] = c;
    }
    c = fgetc( fasta );
  }
  frag_seq->seq[i] = '\0';

  frag_seq->seq_len = i;

  if ( c == '>' ) {
    ungetc( '>', fasta );
    return 1;
  }

  /* Run up against the sequence length limit so truncate it here,
     wind through the fasta filehandle, and return this guy */
  if ( i == INIT_ALN_SEQ_LEN ) {
    while ( (c != '>') &&
	    (c != EOF) ) {
      c = fgetc( fasta );
    }
    if ( c == '>' ) {
      ungetc( '>', fasta );
    }
    fprintf( stderr, "%s is longer than allowed length: %d\n",
	     frag_seq->id, INIT_ALN_SEQ_LEN );
    return 1;
  }

  return 1;
}



/* Return 1 success
 0 failure
 */
int read_fasta_ref(RefSeqP ref, const char* fn) {
  int head_done = 0;
  char c;
  int len, i;
  FILE* ref_f;

  ref->seq = (char*)save_malloc(INIT_REF_SEQ_LEN*sizeof(char));
  ref->size = INIT_REF_SEQ_LEN;
  if (ref->seq == NULL) {
    return 0;
  }

  ref_f = fileOpen(fn, "r");
  if (ref_f == NULL)
      return 0;


  c = fgetc(ref_f);
  if (c == EOF)
    return 0;
  if (c != '>')
    return 0;

  // get id
  len = 0;
  while ( (!isspace( c=fgetc( ref_f ) )) && !head_done) {
    if (c == EOF) {
      return 0;
    }
    ref->id[len] = c;
    len++;
    if (len == MAX_ID_LEN) {
      //id is too long - truncate it
      ref->id[len] = '\0';
      head_done = 1;
    }
  }
  ref->id[len] = '\0';

  // Reset len and head_done to get description
  len = 0;
  head_done = 0;

  /* Check if there is no description */
  if (c == '\n') {
    head_done = 1;
  }
  // If not, don't care about that crappy whitespace we
  // just got
  else {
    c = fgetc(ref_f);
  }
  // But everything else is description
  while ( (c != '\n') && !head_done) {
    if (c == EOF) {
      return 0;
    }
    ref->desc[len] = c;
    len++;
    if (len == MAX_DESC_LEN) {
      // desc is too long, truncate it
      ref->desc[len] = '\0';
      head_done = 1;
    }
    c = fgetc(ref_f);
  }
  ref->desc[len] = '\0';

  // read sequence
  len = 0;
  c = fgetc(ref_f);
  while ( (c != '>' ) && (c != EOF )) {
    if ( isspace(c)) {
      c = fgetc(ref_f);
    }
    else {
      ref->seq[len] = c;
      len++;

      if ( !(len < ref->size)) {
	ref->seq = grow_seq(ref->seq, ref->size);
	ref->size = ref->size * 2;
      }
      c = fgetc(ref_f);
    }
  }
  fclose(ref_f);

  ref->seq[ len ] = '\0';

  /* Now, make and initialize ref->gaps for remembering where
     and how many gaps in the alignment  - NOPE, NOW DONE
     ELSEWHERE - after we learn how long the WRAPPED (!) 
     SEQUENCE IS !!!
  ref->gaps = (int*)save_malloc((len+1)*sizeof(int));
  for (i = 0; i < len; i++) {
    ref->gaps[i] = 0;
    } */

  ref->seq_len = len;

  /* Now, make reverse complement */
  ref->rcseq = (char*)save_malloc(ref->size * sizeof(char));
  if (ref->rcseq == NULL) {
    fprintf( stderr, "Not enough memories for revcom of reference\n");
    exit( 1);
  }
  for (i = 0; i < ref->seq_len; i++) {
    ref->rcseq[i] = revcom_char(ref->seq[ref->seq_len - (i+1)]);
  }
  ref->rcseq[ref->seq_len] = '\0';
  return 1;
}


/* Reads in a set of scoring matrices for each of the PSSM_DEPTH
   positions at the beginning and end of the sequence that are
   to have special scoring matrices and the single 'MIDDLE' matrix
   for everything in the middle.
   Puts these matrices into a PSSM structure.
   Returns a pointer to this structure (PSSMP)
*/
PSSMP read_pssm( const char* fn ) {
  FILE* MF;
  PSSMP submat;
  char* line;
  int cur_pos = 0;
  int i, base;

  line = (char*)save_malloc((MAX_LINE_LEN+1) * sizeof(char));
  MF = fileOpen( fn, "r" );

  if ( MF == NULL ) {
    fprintf( stderr, "Sadly, the substitution matrix cannot be read. Bye bye.\n" );
    exit( 10 );
  }

  /* Initialize the PSSM */
  submat = (PSSMP)save_malloc(sizeof(PSSM));
  submat->depth = PSSM_DEPTH;

  /* Read in the beginning matrices */
  for( cur_pos = 0; cur_pos < PSSM_DEPTH; cur_pos++ ) {
    fgets( line, MAX_LINE_LEN, MF );
    if ( strstr( line, "# Matrix for position" ) == NULL ) {
      fprintf( stderr, "Problem parsing matrix file: %s\n",
	       fn );
      exit( 2 );
    }

    for ( base = 0; base <=3; base++ ) {
      fgets( line, MAX_LINE_LEN, MF );
      sscanf( line, "%d\t%d\t%d\t%d",
	      &submat->sm[cur_pos][base][0],
	      &submat->sm[cur_pos][base][1],
	      &submat->sm[cur_pos][base][2],
	      &submat->sm[cur_pos][base][3] );
      submat->sm[cur_pos][base][4] = N_SCORE; // score for any non-base
    }
    /* ROW of scores for non-base in ref sequence */
    for ( base = 0; base <= 4; base++ ) {
      submat->sm[cur_pos][4][base] = NR_SCORE;
    }

    fgets( line, MAX_LINE_LEN, MF ); // skip blank line
  }

  /* Read in the MIDDLE matrix */
  fgets( line, MAX_LINE_LEN, MF );
  if ( strstr( line, "# Matrix for position: MIDDLE" ) == NULL ) {
    fprintf( stderr, "Problem parsing matrix file, MIDDLE is not where expected: %s\n",
	     fn );
    exit( 2 );
  }
  for ( base = 0; base <=3; base++ ) {
    fgets( line, MAX_LINE_LEN, MF );
    sscanf( line, "%d\t%d\t%d\t%d",
	    &submat->sm[cur_pos][base][0],
	    &submat->sm[cur_pos][base][1],
	    &submat->sm[cur_pos][base][2],
	    &submat->sm[cur_pos][base][3] );
    submat->sm[cur_pos][base][4] = N_SCORE; // 0 score for any non-base
  }
  /* ROW of scores for non-base in ref sequence */
  for ( base = 0; base <= 4; base++ ) {
    submat->sm[cur_pos][4][base] = NR_SCORE;
  }
  fgets( line, MAX_LINE_LEN, MF ); // skip blank line

  /* Read in the ending matrices */
  for( cur_pos = PSSM_DEPTH+1; cur_pos <= (2*PSSM_DEPTH); cur_pos++ ) {
    fgets( line, MAX_LINE_LEN, MF );
    if ( strstr( line, "# Matrix for position:" ) == NULL ) {
      fprintf( stderr, "Problem parsing matrix file: %s\n",
	       fn );
      exit( 2 );
    }
    for ( base = 0; base <=3; base++ ) {
      fgets( line, MAX_LINE_LEN, MF );
      sscanf( line, "%d\t%d\t%d\t%d",
	      &submat->sm[cur_pos][base][0],
	      &submat->sm[cur_pos][base][1],
	      &submat->sm[cur_pos][base][2],
	      &submat->sm[cur_pos][base][3] );
      submat->sm[cur_pos][base][4] = N_SCORE; // 0 score for any non-base
    }
    /* ROW of scores for non-base in ref sequence */
    for ( base = 0; base <= 4; base++ ) {
      submat->sm[cur_pos][4][base] = NR_SCORE;
    }
    fgets( line, MAX_LINE_LEN, MF ); // skip blank line
  }

  free( line );
  fclose( MF );
  submat->depth = PSSM_DEPTH;
  return submat;
}

//TODO: Is this still used?
/* Reads one pairwise alignment from an Udo Stenzel align
 output file of semi-global alignments against a common
 target sequence (usually chrM) into a PWAlnFrag.
 Args: FILE* advanced to next pairwise alignment
 PWAlnFragP to be populated
 Returns 1 if success;
 0 if EOF or failure
 -1 for failure
 */
int read_align_aln(FILE* align_f, PWAlnFragP af) {
  char c;
  int len, strand, aln_len, i;
  int start_gaps = 0;
  int end_gaps = 0;
  int head_done = 0;

  /* Skip past anything until we see a line that begins
     with a '>' */
  c = fgetc(align_f);
  while (c != '>') {
    if (c == EOF) {
      return 0;
    }
    while ( (c != '\n') && (c != EOF)) { // Skip this non >-beginning line
      c = fgetc(align_f);
    }
    c = fgetc(align_f);
  }

  // get ref_id
  len = 0;
  while ( (!isspace( c=fgetc( align_f ) )) && !head_done) {
    if (c == EOF) {
      return 0;
    }
    af->ref_id[len] = c;
    len++;
    if (len == MAX_ID_LEN) {
      //id is too long - truncate it
      af->ref_id[len] = '\0';
      head_done = 1;
    }
  }
  af->ref_id[len] = '\0';

  // Reset len and head_done to get description
  len = 0;
  head_done = 0;

  // Don't care about that crappy whitespace we
  // just got
  c = fgetc(align_f);
  // But everything else is description
  while ( (c != '\n') && !head_done) {
    if (c == EOF) {
      return 0;
    }
    af->ref_desc[len] = c;
    len++;
    if (len == MAX_DESC_LEN) {
      // desc is too long, truncate it
      af->ref_desc[len] = '\0';
      head_done = 1;
    }
    c = fgetc(align_f);
  }
  af->ref_desc[len] = '\0';

  // read ref_seq
  len = 0;
  c = fgetc(align_f);
  while ( (c != '>' ) && (c != EOF )) {
    if (c == '\n' || c == ' ') {
      c = fgetc(align_f);
    } else {
      c = toupper(c);
      af->ref_seq[len] = c;
      len++;

      if (len > INIT_ALN_SEQ_LEN) {
	fprintf( stderr, "Aligned sequence %s is too big\n",
		 af->ref_id);
	return 0; // Too freakin big
      }
      c = fgetc(align_f);
    }
  }
  af->ref_seq[ len ] = '\0';

  if (c == '>')
    ungetc( '>', align_f);

  /* Now, get aligned fragment (frag_id, frag_desc, frag_seq) */
  /* Skip past anything until we see a line that begins
     with a '>' */
  c = fgetc(align_f);
  if (c != '>') {
    if (c == EOF) {
      return 0;
    }
    while ( (c != '\n') && (c != EOF)) { // Skip this non >-beginning line
      c = fgetc(align_f);
    }
    c = fgetc(align_f);
  }

  // get frag_id
  len = 0;
  while ( (!isspace( c=fgetc( align_f ) )) && !head_done) {
    if (c == EOF) {
      return 0;
    }
    af->frag_id[len] = c;
    len++;
    if (len == MAX_ID_LEN) {
      //id is too long - truncate it
      af->frag_id[len] = '\0';
      head_done = 1;
    }
  }
  af->frag_id[len] = '\0';

  // Reset len and head_done to get description
  len = 0;
  head_done = 0;

  // Don't care about that crappy whitespace we
  // just got
  c = fgetc(align_f);
  // But everything else is description
  while ( (c != '\n') && !head_done) {
    if (c == EOF) {
      return 0;
    }
    af->frag_desc[len] = c;
    len++;
    if (len == MAX_DESC_LEN) {
      // desc is too long, truncate it
      af->frag_desc[len] = '\0';
      head_done = 1;
    }
    c = fgetc(align_f);
  }
  af->frag_desc[len] = '\0';

  // read frag_seq
  len = 0;
  c = fgetc(align_f);
  // how much gapped sequence at beginning (context sequence)?
  while (c == '-') {
    start_gaps++;
    af->frag_seq[len++] = c;
    c = fgetc(align_f);
    if (c == '\n' || c == ' ') {
      c = fgetc(align_f);
    }
  }
  while ( (c != '>' ) && (c != EOF )) {
    if (c == '\n' || c == ' ') {
      c = fgetc(align_f);
    } else {
      c = toupper(c);
      af->frag_seq[len++] = c;

      if (len > INIT_ALN_SEQ_LEN) {
	fprintf( stderr, "Aligned sequence %s is too big\n",
		 af->frag_id);
	return 0; // Too freakin big
      }
      // how much gapped sequence at ending (context sequence)?
      if (c == '-') {
	end_gaps++;
      } else {
	end_gaps = 0;
      }
      c = fgetc(align_f);
    }
  }
  af->frag_seq[ len ] = '\0';

  if (c == '>')
    ungetc( '>', align_f);

  // Check length of aligned sequence
  if ( !(strlen(af->frag_seq) == strlen(af->ref_seq) )) {
    fprintf( stderr, "Cannot use %s: ref and frag alignments are unequal lengths\n",
	     af->frag_id);
    // Set score negative so it won't be used
    af->score = -1;
    return 1;
  }

  // Get start, end, strand from pwaln->ref_desc
  if ( !ses_from_align_desc(af, &strand) ) {
    fprintf( stderr, "Problem getting start, end, strand from %s\n",
	     af->ref_desc);
    exit( 1);
  }

  /* Learn whether adapter was cut off or not from af->frag_desc
     and set af->trimmed accordingly */
  if ( !adapt_from_desc(af) ) {
    fprintf( stderr, "Problem learning from %s if adapter was cut, set to not\n",
	     af->frag_desc);
  }

  // Do reverse complement of both if minus strand alignment
  if (strand == -1) {
    revcom_PWAF(af);
    af->revcom = 1;
  } else {
    af->revcom = 0;
  }

  /* Now, remove context from reference, if any, and update
     start and end coordinates */
  aln_len = strlen(af->ref_seq);
  aln_len -= start_gaps;
  aln_len -= end_gaps;

  if (af->revcom) {
    for (i = 0; i < aln_len; i++) {
      af->ref_seq[i] = af->ref_seq[i+end_gaps];
      af->frag_seq[i] = af->frag_seq[i+end_gaps];
    }
    af->start += end_gaps;
    af->end -= start_gaps;
  } else {
    for (i = 0; i < aln_len; i++) {
      af->ref_seq[i] = af->ref_seq[i+start_gaps];
      af->frag_seq[i] = af->frag_seq[i+start_gaps];
    }
    af->start += start_gaps;
    af->end -= end_gaps;
  }
  af->ref_seq[i] = '\0';
  af->frag_seq[i] = '\0';

  /* check if the reference sequence was complemented; if so, invert the
   * revcom flag.  everything else isn't affected. */
  if (af->frag_desc[0] == '-') {
    af->revcom = !af->revcom;
  }
  return 1;
}





void ace_output(MapAlignmentP maln) {
	int number_of_contigs = 1;
	char* consensus = get_consensus(maln);
	int number_of_BS = 1;
	int number_of_reads = maln->num_aln_seqs;
	const int QUALITY_SCORE = 40;
	int number_bases = get_consensus_length(maln);
	char* contig_name = maln->ref->id;

	int i, j, line_pos = 0, gaps;

	int max_line_length = 50;
	AlnSeqP aln_seq;

	//////////////////////////////////////////// ASSEMBLY INFORMATION (AS) ////////////////////////////////////////////////
	printf("AS %d %d\n\n", number_of_contigs, number_of_reads + 1); // if we allow repairing we have one (fake) read more

	//////////////////////////////////////////// CONTIG INFORMATION (CO) //////////////////////////////////////////////////
	printf("CO %s %d %d %d %c\n", contig_name, number_bases, number_of_reads + 1,
			number_of_BS, 'U');

	//////////////////////////////////////////// CONSENSUS  ///////////////////////////////////////////////////////////////
	// print consenus --> padded positions must be * in ace

	char curr_line[max_line_length + 1];
	for (i = 0; i < number_bases; i++) {

		curr_line[line_pos++] = (consensus[i] == '-') ? '*'
					: ( (consensus[i] == ' ') ? ('X') : (consensus[i]) );

		if (line_pos == max_line_length) {
			curr_line[line_pos] = '\0';
			printf("%s\n", curr_line);
			line_pos = 0;
		}
	}
	curr_line[line_pos] = '\0';
	printf("%s\n", curr_line);
	printf("\n");

	//////////////////////////////////////////////// BASE QUALITIES //////////////////////////////////////////
	// Print quality scores (UNPADDED!)
	printf("BQ\n");
	for (i = 0; i < number_bases; i++) {
		if (consensus[i] != '-')
			printf("%d ", QUALITY_SCORE); // mia has no quality scores -> Maybe I can build something from the coverage informations ...
		if (i % max_line_length == 0)
			printf("\n");
	}
	printf("\n\n");

	/////////////////////////////////////////////// ORDER OF THE READS (AF) (PADDED!) ////////////////////////////

	printf("AF FAKE_READ-IGNORE_ME U %d\n", 1);
	for (i = 0; i < number_of_reads; i++) {
		aln_seq = maln->AlnSeqArray[i];
		gaps = sum_of_gaps(maln, aln_seq->start); // How many gaps upstream?
		printf("AF %s %c %d\n", aln_seq->id, (aln_seq->revcom) ? 'C' : 'U',
				aln_seq->start + gaps+1);
	}

	printf("\n");


	////////////////////////////////////////////// BASE SEGMENTS (BS) ////////////////////////////////////////////
	printf("BS 1 %d %s\n", strlen(consensus), "FAKE_READ-IGNORE_ME");
	printf("\n");

	///////////////////////////////////////////////// PRINT THE READS  (RD) ////////////////////////////////////////
	int tmp;
        maln->ref->gaps[maln->ref->seq_len] = 0;
	for (tmp = 0; tmp < number_of_reads; tmp++) {
		aln_seq = maln->AlnSeqArray[tmp];
		int gaps = 0;
                int n_gaps = 0;
                int ix = 0;
                int ins_len = 0;
                char * ins = NULL;

                //if (aln_seq->end >= maln->ref->seq_len)
                //    aln_seq->end = maln->ref->seq_len - 1;

               // printf("%d %d %d\n", aln_seq->end, maln->ref->seq_len , maln->ref->gaps[aln_seq->end]);

                /* Find how many gaps are in this region */
                for (i = aln_seq->start; i <= aln_seq->end; i++) {
                    gaps += maln->ref->gaps[i];
                }

		printf("RD %s %d %d %d\n", aln_seq->id,
				strlen(aln_seq->seq) + gaps, 0, 0);
		char *sequence = (char*)save_malloc(strlen(aln_seq->seq) + gaps + 1
				* sizeof(char));
		j=0;
		for (i = aln_seq->start; i<= aln_seq->end; i++) {
			if (maln->ref->gaps[i] > 0) {
                            if (aln_seq->ins[i - aln_seq->start] != NULL){
                                ins = aln_seq->ins[i - aln_seq->start];
                                ins_len = strlen(ins);
                            }
                            else
                                ins_len = 0;
                            for (n_gaps = 0; n_gaps < maln->ref->gaps[i]; n_gaps++) {
                                if (n_gaps < ins_len)
                                    sequence[j++] = ins[n_gaps];
                                else
                                    sequence[j++] = '*';
				}

			}
			sequence[j++] = aln_seq->seq[i - aln_seq->start];
		}

		line_pos = 0;
		for (i = 0; i < j; i++) {
			curr_line[line_pos++] = (sequence[i] == '-')?'*':sequence[i];

			if (line_pos == max_line_length) {
				curr_line[line_pos] = '\0';
				printf("%s\n", curr_line);
				line_pos = 0;
			}
		}
		curr_line[line_pos] = '\0';
		printf("%s\n", curr_line);
		printf("\n");

		printf("QA %d %d %d %d\n", 1, strlen(aln_seq->seq) + gaps, 1,
				strlen(aln_seq->seq) + gaps);

		printf(
				"DS CHROMAT_FILE: %s PHD_FILE: %s_FAKE.phd TIME: Tue Feb 21 15:42:35 1984\n\n",
				aln_seq->id, aln_seq->id);

	}

		printf("RD FAKE_READ-IGNORE_ME %d %d %d\n", number_bases, 0, 0);

		line_pos= 0;
			for (i = 0; i < number_bases; i++) {

				curr_line[line_pos++] = (consensus[i] == '-') ? '*'
							: ( (consensus[i] == ' ') ? ('X') : (consensus[i]) );

				if (line_pos == max_line_length) {
					curr_line[line_pos] = '\0';
					printf("%s\n", curr_line);
					line_pos = 0;
				}
			}
			curr_line[line_pos] = '\0';
			printf("%s\n", curr_line);
		printf("\n\n");
		printf("QA %d %d %d %d\n", 1, number_bases, 1, number_bases);
		printf("DS CHROMAT_FILE: %s PHD_FILE: %s_FAKE.phd TIME: Tue Feb 21 23:23:23 1984\n", "FAKE_READ", "FAKE_READ");


	return;
}



/** fileOpen **/
FILE * fileOpen(const char *name, char access_mode[]) {
	FILE * f;
	f = fopen(name, access_mode);
	if (f == NULL) {
		fprintf( stderr, "%s\n", name);
		perror("Cannot open file");
		return NULL;
	}
	return f;
}

void fasta_print_cons(char* cons, char* id) {
	int len, i, line_pos;
	char curr_line[FASTA_LINE_WIDTH + 1];
	len = strlen(cons);
	printf(">%s\n", id);
	line_pos = 0;
	for (i = 0; i < len; i++) {
		if ( !(cons[i] == '-')) {
			if (cons[i] == ' ') {
				curr_line[line_pos++] = 'X';
			} else {
				curr_line[line_pos++] = cons[i];
			}
			if (line_pos == FASTA_LINE_WIDTH) {
				curr_line[line_pos] = '\0';
				printf("%s\n", curr_line);
				line_pos = 0;
			}
		}
	}
	curr_line[line_pos] = '\0';
	printf("%s\n", curr_line);
}

void fasta_aln_print(char* seq, char* id) {
	int len, i, line_pos;
	char curr_line[FASTA_LINE_WIDTH + 1];
	len = strlen(seq);
	printf(">%s\n", id);
	line_pos = 0;
	for (i = 0; i < len; i++) {
		if (seq[i] == ' ') {
			curr_line[line_pos++] = 'X';
		} else {
			curr_line[line_pos++] = seq[i];
		}
		if (line_pos == FASTA_LINE_WIDTH) {
			curr_line[line_pos] = '\0';
			printf("%s\n", curr_line);
			line_pos = 0;
		}
	}
	curr_line[line_pos] = '\0';
	printf("%s\n", curr_line);
}


void clustalw_print_cons(char* cons, char* aln_ref, char* ref_id) {
	int len, ln, i, ln_len, ref_id_len;
	char ref_start[18];
	char* curr_ref_line;
	char* curr_cons_line;
	curr_ref_line = (char*)save_malloc((CLUSTALW_LINE_WIDTH+1) * sizeof(char));
	curr_cons_line = (char*)save_malloc((CLUSTALW_LINE_WIDTH+1) * sizeof(char));
	len = strlen(cons);
	ref_id_len = strlen(ref_id);
	ln = 0;

	strncpy(ref_start, ref_id, 15);

	for (i = ref_id_len; i < 15; i++) {
		ref_start[i] = ' ';
	}
	ref_start[15] = ' ';
	ref_start[16] = ' ';
	ref_start[17] = '\0';

	printf("CLUSTAL W (1.8) multiple sequence alignment\n");
	while ( (ln * CLUSTALW_LINE_WIDTH) < len) {
		// First make and print the Reference sequence line
		strncpy(curr_ref_line, &aln_ref[CLUSTALW_LINE_WIDTH*ln], 60);
		curr_ref_line[CLUSTALW_LINE_WIDTH] = '\0';
		printf("%s%s\n", ref_start, curr_ref_line);

		// Then, the consensus
		strncpy(curr_cons_line, &cons[CLUSTALW_LINE_WIDTH*ln], 60);
		curr_cons_line[CLUSTALW_LINE_WIDTH] = '\0';
		// Replace spaces with X, to denote no coverage
		for (i = 0; i < CLUSTALW_LINE_WIDTH; i++) {
			if (curr_cons_line[i] == ' ') {
				curr_cons_line[i] = 'X';
			}
		}
		printf("Consensus        %s\n", curr_cons_line);

		// Then, the * line where there is agreement
		ln_len = strlen(curr_cons_line);
		printf("                 ");
		for (i = 0; i < ln_len; i++) {
			if (curr_ref_line[i] == curr_cons_line[i]) {
				printf("*");
			} else {
				printf(" ");
			}
		}
		printf("\n\n\n");
		ln++;
	}
	free(curr_ref_line);
	free(curr_cons_line);
}


void line_print_cons(char* consensus, char* aln_ref, char* ref_id, int* cov) {
	int len, i;
	len = strlen(consensus);
	printf("Consensus, %s, coverage:\n", ref_id);
	printf("%s\n%s\n", consensus, aln_ref);

	for (i = 0; i < len; i++) {
		printf("%d ", cov[i]);
	}
	printf("\n");
}

void color_print(char* string) {
	char c = *string++;
	while (c != '\0') {
		switch (c) { // insert color codes. see: http://linuxgazette.net/issue65/padala.html
		// Green background = 42
		case 'a':
			;
		case 'A':
			printf("\33[37;42m");
			break; //4
			// Blue background = 44
		case 'c':
			;
		case 'C':
			printf("\33[37;44m");
			break; //4
			// Black background = 40
		case 'g':
			;
		case 'G':
			printf("\33[37;40m");
			break; //4
			// Red background = 41
		case 't':
			;
		case 'T':
			printf("\33[37;41m");
			break; //4
			// Gray background = 47
		case '-':
			printf("\33[47;30m");
			break; //4
			// white = 47
		default:
			printf("\33[0m");
		}

		printf("%c", c);
		c = (*string++);
	}
	printf("\33[0m\n");
}


IDsListP parse_ids(char* fn) {
	int i;
	int id_num = 0;
	int c_num = 0;
	IDsListP ids;
	FILE* IDS;
	char** ids_array;
	char* first_id;
	char c;
	// First, allocate the IDsList
	ids = (IDsListP)save_malloc(sizeof(IDsList));

	ids_array = (char**)save_malloc(INIT_NUM_IDS * sizeof( char* ));
	first_id = (char*)save_malloc(INIT_NUM_IDS * MAX_ID_LEN * sizeof(char));

	for (i = 0; i < INIT_NUM_IDS; i++) {
		ids_array[i] = &first_id[i*MAX_ID_LEN];
	}
	ids->ids = ids_array;
	ids->size = INIT_NUM_IDS;
	ids->sorted = 0;
	/* Parse the file, loading up ids_array with IDs */
	IDS = fileOpen(fn, "r");
	c = fgetc(IDS);
	while (c != EOF) {
		if (c == '\n') { // just finished a line (ID)
			ids->ids[id_num++][c_num] = '\0';
			c_num = 0;
			if (id_num >= ids->size) {
				grow_ids_list(ids);
			}
		} else {
			if (c_num >= MAX_ID_LEN) {
				ids->ids[id_num][c_num] = '\0';
			} else {
				ids->ids[id_num][c_num++] = c;
			}
		}
		c = fgetc(IDS);
	}
	fclose(IDS);

	ids->num_ids = id_num;
	qsort(ids->ids, id_num, sizeof(char* ), idCmp);
	ids->sorted = 1;

	return ids;
}


