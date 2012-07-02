#include "map_alignment.h"

/* Initialize a MapAlignment object and return a pointer to it */
MapAlignmentP init_map_alignment(void) {
    MapAlignmentP aln;
    AlnSeqP first_seq, as;
    size_t i, j;

    // First, allocate the alignment
    aln = (MapAlignmentP) save_malloc(sizeof (MapAlignment));

    // Allocate memory for the RefSeq
    aln->ref = (RefSeqP) save_malloc(sizeof (RefSeq));
    if (aln->ref == NULL) {
        return NULL;
    }
    // Zero-out the RefSeq
    for (i = 0; i <= MAX_ID_LEN; i++) {
        aln->ref->id[i] = '\0';
    }
    for (i = 0; i <= MAX_DESC_LEN; i++) {
        aln->ref->desc[i] = '\0';
    }
    aln->ref->seq = NULL;
    aln->ref->rcseq = NULL;
    aln->ref->size = 0;
    aln->ref->gaps = NULL;
    aln->ref->circular = 0;
    aln->ref->wrap_seq_len = 0;

    // Allocate memory for all the aligned sequences
    first_seq = (AlnSeqP) save_malloc(INIT_NUM_ALN_SEQS *
            sizeof ( AlnSeq));
    if (first_seq == NULL) {
        return NULL;
    }

    // Now, allocate the array of pointers to the
    // aligned seqs
    aln->AlnSeqArray = (AlnSeqP*) save_malloc(INIT_NUM_ALN_SEQS *
            sizeof ( AlnSeqP));
    if (aln->AlnSeqArray == NULL) {
        return NULL;
    }

    // Now, point the pointers to the pointees
    for (i = 0; i < INIT_NUM_ALN_SEQS; i++) {
        aln->AlnSeqArray[i] = &first_seq[i];
        /* Zero them out */
        as = aln->AlnSeqArray[i];
        for (j = 0; j <= MAX_ID_LEN; j++) {
            as->id[j] = '\0';
        }
        for (j = 0; j <= MAX_DESC_LEN; j++) {
            as->desc[j] = '\0';
        }
        for (j = 0; j <= INIT_ALN_SEQ_LEN; j++) {
            as->seq[j] = '\0';
        }
        /* Set all their char* ins to NULL */
        for (j = 0; j <= INIT_ALN_SEQ_LEN; j++) {
            as->ins[j] = NULL;
        }
        as->start = 0;
        as->end = 0;
        as->revcom = 0;
        as->trimmed = 0;
        as->score = 0;
        as->segment = 'n';
    }

    aln->size = INIT_NUM_ALN_SEQS;
    aln->num_aln_seqs = 0;

    return aln;
}

/* free_map_alignment
 Takes a MapAlignmentP (maln)
 Frees the memory pointed to by its components
 Returns nothing
 */
void free_map_alignment(MapAlignmentP maln) {
    /* First, free the maln->ref components */
    free(maln->ref->seq);
    if (maln->ref->rcseq != NULL) {
        free(maln->ref->rcseq);
    }
    free(maln->ref->gaps);
    free(maln->ref);

    /* Now, free the AlnSeqArray */
    free(maln->AlnSeqArray[0]);
    free(maln->AlnSeqArray);

    /* Now, free the MapAlignment */
    free(maln);

    /* That oughta do it */
    return;
}

void show_consensus(MapAlignmentP maln, int out_format) {
    char* consensus;
    char* aln_ref;
    char* ins_cons;
    char cons_id[MAX_ID_LEN + 1];
    int i, j, cons_pos, ref_pos, ref_gaps;
    int* cov;
    int* ins_cov;
    int* ref_poss;
    int len_consensus = get_consensus_length(maln);
    AlnSeqP aln_seq;
    BaseCountsP bcs;
    PSSMP psm;

    bcs = (BaseCountsP) save_malloc(sizeof (BaseCounts));
    reset_base_counts(bcs);

    consensus = (char*) save_malloc((len_consensus + 1) * sizeof (char));
    aln_ref = (char*) save_malloc((len_consensus + 1) * sizeof (char));
    cov = (int*) save_malloc((len_consensus + 1) * sizeof (int));
    ref_poss = (int*) save_malloc((len_consensus + 1) * sizeof (int));

    ins_cons = (char*) save_malloc(MAX_INS_LEN * sizeof (char));
    ins_cov = (int*) save_malloc(MAX_INS_LEN * sizeof (int));

    cons_pos = 0;
    ref_pos = 0;
    /* Go through each position of the reference sequence */
    for (ref_pos = 0; ref_pos < maln->ref->seq_len; ref_pos++) {
        /* How many gaps preceeded this position? */
        ref_gaps = maln->ref->gaps[ref_pos];

        /* Add these gaps to the reference aligned string */
        if ((ref_gaps > 0) && (ref_pos > 0)) {
            find_ins_cons(maln, ref_pos, ins_cons, ins_cov, out_format);
            for (j = 0; j < ref_gaps; j++) {
                aln_ref[cons_pos] = '-';
                consensus[cons_pos] = ins_cons[j];
                cov[cons_pos] = ins_cov[j];
                ref_poss[cons_pos] = ref_pos;
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
            if ((aln_seq->start <= ref_pos) && // checked
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
        cov[cons_pos] = bcs->cov;
        ref_poss[cons_pos] = ref_pos;
        if ((out_format == 4) && !(aln_ref[cons_pos] == consensus[cons_pos])) {
            show_single_pos(ref_pos, aln_ref[cons_pos], consensus[cons_pos],
                    bcs);
        }
        if (out_format == 41) {
            show_single_pos(ref_pos, aln_ref[cons_pos], consensus[cons_pos],
                    bcs);
        }
        cons_pos++;
    }
    consensus[cons_pos] = '\0';
    aln_ref[cons_pos] = '\0';

    /* Now, output the reference and consensus sequences and the
       coverage in specified way */
    switch (out_format) {
        case 1:
            clustalw_print_cons(consensus, aln_ref, maln->ref->id);
            break;
        case 2:
            line_print_cons(consensus, aln_ref, maln->ref->id, cov);
            break;
        case 3:
            /* Add starts and ends info */
            print_assembly_summary(maln);
            col_print_cons(consensus, aln_ref, cov, ref_poss, maln);
            break;
        case 4:
            ; /* Do nothing, this one is checked along the way */
            break;
        case 41:
            ; /* Do nothing, this one is checked along the way */
            break;
        case 5:
            //		sprintf(cons_id, "%s-assembled", maln->ref->id);
            fasta_print_cons(consensus, maln->ref->id);
            break;
    }

    /* Free memory! */
    free(bcs);
    free(consensus);
    free(aln_ref);
    free(cov);
    free(ref_poss);
    free(ins_cons);
    free(ins_cov);
}

int get_consensus_length(MapAlignmentP maln) {
    int i, num_gaps = 0;
    for (i = 0; i < maln->ref->seq_len; i++)
        num_gaps += maln->ref->gaps[i];
    return maln->ref->seq_len + num_gaps;
}

char* get_consensus(MapAlignmentP maln) {
    int len_consensus = get_consensus_length(maln);
    char* consensus = (char*) save_malloc((len_consensus + 1) * sizeof (char));
    char* ins_cons = (char*) save_malloc(MAX_INS_LEN * sizeof (char));
    int i, j, cons_pos, ref_pos, ref_gaps;    
    int* ins_cov = (int*) save_malloc(MAX_INS_LEN * sizeof (int));    
    AlnSeqP aln_seq;
    BaseCountsP bcs;
    PSSMP psm;
    bcs = (BaseCountsP) save_malloc(sizeof (BaseCounts));
    reset_base_counts(bcs);

    cons_pos = 0;
    ref_pos = 0;
    /* Go through each position of the reference sequence */
    for (ref_pos = 0; ref_pos < maln->ref->seq_len; ref_pos++) {
        /* How many gaps preceeded this position? */
        ref_gaps = maln->ref->gaps[ref_pos];

        /* Add these gaps to the reference aligned string */
        if ((ref_gaps > 0) && (ref_pos > 0)) {
            find_ins_cons(maln, ref_pos, ins_cons, ins_cov, 5);
            for (j = 0; j < ref_gaps; j++) {
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
            if ((aln_seq->start <= ref_pos) && // checked
                    (aln_seq->end >= ref_pos)) {

                psm = (aln_seq->revcom) ? (maln->rpsm) : (maln->fpsm);

                add_base(aln_seq->seq[ref_pos - aln_seq->start], bcs, psm,
                        aln_seq->smp[ref_pos - aln_seq->start]);
            }
        }
        consensus[cons_pos] = find_consensus(bcs, maln->cons_code);
        cons_pos++;
    }
    consensus[cons_pos] = '\0';
    return consensus;
}

/* Write out the data in a MapAlignment data structure
 to a file
 */
int write_ma(char* fn, MapAlignmentP maln) {
    int i, j, row, col;
    //    char* at;
    time_t t;
    int aln_seq_len;
    FILE* MAF;
    AlnSeqP as;
    PSSMP fpsm, rpsm;

    MAF = fileOpen(fn, "w");

    t = time(NULL);
    //at = (char*) save_malloc(64 * sizeof (char));

    //at = asctime(localtime(&t));
    /* First, write a nice header */
    fprintf(MAF, "/* map_alignment [V%s] */ %s",PACKAGE_VERSION , 
	    asctime(localtime(&t)) );

    /* Write MapAlignment Info */
    fprintf(MAF, "MALN_NAS %d\n", maln->num_aln_seqs);
    fprintf(MAF, "MALN_SIZ %d\n", maln->size);
    fprintf(MAF, "MALN_COC %d\n", maln->cons_code);

    /* Write the reference sequence and associated data */
    fprintf(MAF, "__REFERENCE__\n");
    fprintf(MAF, "ID %s\n", maln->ref->id);
    fprintf(MAF, "DESC %s\n", maln->ref->desc);
    fprintf(MAF, "LEN %d\n", maln->ref->seq_len);
    fprintf(MAF, "SIZE %d\n", maln->ref->size);
    fprintf(MAF, "SEQ ");
    for (i = 0; i < maln->ref->seq_len; i++) {
        fprintf(MAF, "%c", maln->ref->seq[i]);
    }
    fprintf(MAF, "\n");

    fprintf(MAF, "GAPS");
    for (i = 0; i < maln->ref->seq_len; i++) {
        fprintf(MAF, " %d", maln->ref->gaps[i]);
    }
    fprintf(MAF, "\n");

    /* Write the PSSMs */
    fpsm = maln->fpsm;
    rpsm = maln->rpsm;
    fprintf(MAF, "__PSSM__\n");
    fprintf(MAF, "DEPTH %d\n", fpsm->depth);
    fprintf(MAF, "FPSM:\n");
    for (i = 0; i <= (fpsm->depth * 2); i++) {
      for (row = 0; row <= 4; row++) {
	fprintf(MAF, "%d %d %d %d %d\n", 
		fpsm->sm[i][row][0],
		fpsm->sm[i][row][1], 
		fpsm->sm[i][row][2],
		fpsm->sm[i][row][3], 
		fpsm->sm[i][row][4]);
      }
      fprintf(MAF, "\n");
    }
    
    fprintf(MAF, "RPSM:\n");
    for (i = 0; i <= (fpsm->depth * 2); i++) {
      for (row = 0; row <= 4; row++) {
	fprintf(MAF, "%d %d %d %d %d\n", rpsm->sm[i][row][0],
		rpsm->sm[i][row][1], rpsm->sm[i][row][2],
		rpsm->sm[i][row][3], rpsm->sm[i][row][4]);
      }
      fprintf(MAF, "\n");
    }

    /* Write all the aligned fragments */
    fprintf(MAF, "__ALNSEQS__\n");
    for (i = 0; i < maln->num_aln_seqs; i++) {
      as = maln->AlnSeqArray[i];
      aln_seq_len = strlen(as->seq);
      fprintf(MAF, "ID %s\n", as->id);
      fprintf(MAF, "DESC %s\n", as->desc);
      fprintf(MAF, "SCORE %d\n", as->score);
      fprintf(MAF, "NUM_INPUTS %d\n", as->num_inputs);
      fprintf(MAF, "START %d\n", as->start);
      fprintf(MAF, "END %d\n", as->end);
      fprintf(MAF, "RC %d\n", as->revcom);
      fprintf(MAF, "TR %d\n", as->trimmed);
      fprintf(MAF, "SEG %c\n", as->segment);
      fprintf(MAF, "SEQ %s\n", as->seq);
      fprintf(MAF, "SMP %s\n", as->smp);
      fprintf(MAF, "INS_POS");
        for (j = 0; j < aln_seq_len; j++) {
            if (as->ins[j] == NULL) {
                //	fprintf( MAF, " _" );
            } else {
                fprintf(MAF, " %d %s", j, as->ins[j]);
            }
        }
        fprintf(MAF, "\n");
    }
    fclose(MAF);
    return 1;
}

MapAlignmentP read_ma(const char* fn) {
    MapAlignmentP maln;
    AlnSeqP as;
    FILE* MAF;
    char* line;
    char* tmp_ins;
    char c;
    int tmp, i, as_num, ins_pos, depth, row, A, C, G, T, N;

    line = (char*) save_malloc((MAX_LINE_LEN + 1) * sizeof (char));
    MAF = fileOpen(fn, "r");

    maln = init_map_alignment();
    maln->fpsm = (PSSMP) save_malloc(sizeof (PSSM));
    maln->rpsm = (PSSMP) save_malloc(sizeof (PSSM));

    /* Check header */
    fgets(line, MAX_LINE_LEN, MAF);
    if (strstr(line, "/* map_alignment") == NULL) {
        fprintf(stderr, "%s does not look like a map_alignment input file\n",
                fn);
        exit(1);
    }

    /* Parse MALN_NAS */
    fgets(line, MAX_LINE_LEN, MAF);
    sscanf(line, "MALN_NAS %d", &maln->num_aln_seqs);

    /* Parse MALN_SIZ; grow the AlnSeqArray of the MapAlignment until
     it's at least as big as before */
    fgets(line, MAX_LINE_LEN, MAF);
    sscanf(line, "MALN_SIZ %d", &tmp);
    while (maln->size < tmp) {
        grow_alns_map_alignment(maln);
    }

    /* Parse MALN_NAS */
    fgets(line, MAX_LINE_LEN, MAF);
    sscanf(line, "MALN_COC %d", &maln->cons_code);

    /* Parse the reference sequence header */
    fgets(line, MAX_LINE_LEN, MAF);
    if (strstr(line, "__REFERENCE__") == NULL) {
        fprintf(stderr, "Do not see reference sequence header in %s\n", fn);
        exit(1);
    }

    /* Parse the reference ID */
    fgets(line, MAX_LINE_LEN, MAF);
    sscanf(line, "ID %s", maln->ref->id);

    /* Parse the reference DESC */
    fgets(line, MAX_LINE_LEN, MAF);
    sscanf(line, "DESC %s", maln->ref->desc);

    /* Parse the reference LEN and make the maln->ref->gaps
     point to an array of ints this size*/
    fgets(line, MAX_LINE_LEN, MAF);
    sscanf(line, "LEN %d", &maln->ref->seq_len);
    maln->ref->gaps = (int*) save_malloc(maln->ref->seq_len * sizeof (int));

    /* Parse the reference SIZE and make the maln->ref->seq point
     to a char array this size */
    fgets(line, MAX_LINE_LEN, MAF);
    sscanf(line, "SIZE %d", &maln->ref->size);
    maln->ref->seq = (char*) save_malloc(maln->ref->size * sizeof (char));

    /* Parse the reference SEQ and put it into maln->ref->seq */
    fgets(line, MAX_LINE_LEN, MAF);
    sscanf(line, "SEQ %s", maln->ref->seq);

    /* Check to make sure the LEN info we got a few lines back is correct */
    if (!(strlen(maln->ref->seq) == maln->ref->seq_len)) {
        fprintf(stderr, "Reported length of reference sequence %d is not observed length %d\n",
                (int) maln->ref->seq_len, (int) strlen(maln->ref->seq));
        exit(1);
    }

    /* Parse the reference GAPS and put them into maln->ref->gaps */
    fscanf(MAF, "GAPS"); // Go past the GAPS string and start getting %d
    for (i = 0; i < maln->ref->seq_len; i++) {
        fscanf(MAF, " %u", &maln->ref->gaps[i]);
    }

    /* Yoink-a-doink on past the \n that we never passed from the GAPS line */
    c = fgetc(MAF);
    while ((c != '\n') && (c != EOF)) {
        c = fgetc(MAF);
    }

    /* Parse the line that announces we're starting the PSSM section */
    fgets(line, MAX_LINE_LEN, MAF);
    if (strstr(line, "__PSSM__") == NULL) {
        fprintf(stderr, "Do not see __PSSM__ line in %s\n", fn);
        exit(2);
    }

    /* Get/set the PSSMP depths */
    fgets(line, MAX_LINE_LEN, MAF);
    sscanf(line, "DEPTH %d", &depth);
    maln->fpsm->depth = depth;
    maln->rpsm->depth = depth;

    /* Skip past the FPSM: line */
    fgets(line, MAX_LINE_LEN, MAF);
    if (strstr(line, "FPSM:") == NULL) {
        fprintf(stderr, "Do not see the FPSM: in %s\n", fn);
        exit(2);
    }

    /* Now, get depth number of substitution matrices for the maln->fpsm->sm */
    for (i = 0; i <= (depth * 2); i++) {
        for (row = 0; row <= 4; row++) {
            fgets(line, MAX_LINE_LEN, MAF);
            sscanf(line, "%d %d %d %d %d", &A, &C, &G, &T, &N);
            maln->fpsm->sm[i][row][0] = A;
            maln->fpsm->sm[i][row][1] = C;
            maln->fpsm->sm[i][row][2] = G;
            maln->fpsm->sm[i][row][3] = T;
            maln->fpsm->sm[i][row][4] = N;
        }
        /* Skip blank line separating matrices */
        fgets(line, MAX_LINE_LEN, MAF);
    }

    /* Skip past the RPSM: line */
    fgets(line, MAX_LINE_LEN, MAF);
    if (strstr(line, "RPSM:") == NULL) {
        fprintf(stderr, "Do not see the RPSM: in %s\n", fn);
        exit(2);
    }

    /* Now, get depth number of substitution matrices for the maln->rpsm->sm */
    for (i = 0; i <= (depth * 2); i++) {
        for (row = 0; row <= 4; row++) {
            fgets(line, MAX_LINE_LEN, MAF);
            sscanf(line, "%d %d %d %d %d", &A, &C, &G, &T, &N);
            maln->rpsm->sm[i][row][0] = A;
            maln->rpsm->sm[i][row][1] = C;
            maln->rpsm->sm[i][row][2] = G;
            maln->rpsm->sm[i][row][3] = T;
            maln->rpsm->sm[i][row][4] = N;
        }
        /* Skip blank line separating matrices */
        fgets(line, MAX_LINE_LEN, MAF);
    }

    /* Parse the line that announces we're starting the aligned fragments
     section of output */
    fgets(line, MAX_LINE_LEN, MAF);
    if (strstr(line, "__ALNSEQS__") == NULL) {
        fprintf(stderr, "Do not see __ALNSEQS__ line in %s\n", fn);
        exit(1);
    }

    /* Go through the parsing for as many aligned fragments as we're
     expecting */
    for (as_num = 0; as_num < maln->num_aln_seqs; as_num++) {
        as = maln->AlnSeqArray[as_num];

        /* Get ID line */
        fgets(line, MAX_LINE_LEN, MAF);
        sscanf(line, "ID %s\n", as->id);

        /* Get DESC line */
        fgets(line, MAX_LINE_LEN, MAF);
        strcpy(as->desc, &line[5]);
        as->desc[strlen(as->desc) - 1] = '\0'; // get rid of \n

        /* Get SCORE line */
        fgets(line, MAX_LINE_LEN, MAF);
        sscanf(line, "SCORE %d\n", &as->score);

	/* Get NUM_INPUTS line, if there */
	fgets(line, MAX_LINE_LEN, MAF);
	if ( sscanf( line, "NUM_INPUTS %d\n", &as->num_inputs ) == 1 ) {
	  fgets(line, MAX_LINE_LEN, MAF);
	}
	else {
	  as->num_inputs = 1;
	}

        /* Get START line */
        sscanf(line, "START %d\n", &as->start);

        /* Get END line */
        fgets(line, MAX_LINE_LEN, MAF);
        sscanf(line, "END %d\n", &as->end);

        /* Get RC line */
        fgets(line, MAX_LINE_LEN, MAF);
        sscanf(line, "RC %d\n", &as->revcom);

        /* Get TR line */
        fgets(line, MAX_LINE_LEN, MAF);
        sscanf(line, "TR %d\n", &as->trimmed);

        /* Get SEG line */
        fgets(line, MAX_LINE_LEN, MAF);
        sscanf(line, "SEG %c\n", &as->segment);

        /* Get SEQ line */
        fgets(line, MAX_LINE_LEN, MAF);
        sscanf(line, "SEQ %s\n", as->seq);

        /* Get SMP line */
        fgets(line, MAX_LINE_LEN, MAF);
        sscanf(line, "SMP %s\n", as->smp);

        /* Get INS line */
        tmp_ins = (char*) save_malloc(MAX_INS_LEN * sizeof (char));

        fscanf(MAF, "INS_POS");
        while (fscanf(MAF, " %d %s", &ins_pos, tmp_ins) == 2) {
            as->ins[ins_pos] = (char*) save_malloc(MAX_INS_LEN * sizeof (char));
            strcpy(as->ins[ins_pos], tmp_ins);
        }

    }
    fclose(MAF);
    free(line);
    return maln;
}

int count_aln_seqs(MapAlignmentP maln) {
    int i;
    int tot_aln_seqs = 0;

    /* Count each one that is segment a, n, or f */
    for (i = 0; i < maln->num_aln_seqs; i++) {
        if (maln->AlnSeqArray[i]->segment != 'b') {
            tot_aln_seqs++;
        }
    }
    return tot_aln_seqs;
}

/* Sorts the AlnSeqArray by alnSeqCmp (its start and end
 coordinates.
 Note that after this operation, any FragSeqDB pointing
 to this AlnSeqArray will be wrong! */
void sort_aln_frags(MapAlignmentP maln) {
    qsort((void*) maln->AlnSeqArray, (size_t) maln->num_aln_seqs,
            sizeof (AlnSeqP), alnSeqCmp);
}

void print_assembly_summary(MapAlignmentP maln) {
    int i;
    int total_frag_len = 0;

    for (i = 0; i < maln->num_aln_seqs; i++) {
        total_frag_len += (maln->AlnSeqArray[i]->end
                - maln->AlnSeqArray[i]->start + 1);
    }

    printf("# Map reference ID: %s\n", maln->ref->id);
    printf("# Map reference length: %d\n", maln->ref->seq_len);
    printf("# Number of fragments aligned to reference: %d\n",
            count_aln_seqs(maln));
    //	  maln->num_aln_seqs );
    printf("# Total length of aligned fragments: %d\n", total_frag_len);
    printf("# Average coverage: %0.3f\n", ((double) total_frag_len
            / (double) maln->ref->seq_len));

}


// Return the absolute number of gaps upstream of this position

int sum_of_gaps(MapAlignmentP maln, int pos) {
    int i, gaps;
    gaps = 0;
    for (i = 0; (i < pos); i++)
        gaps += maln->ref->gaps[i];
    return gaps;
}

/* Grow the space for a MapAlignment to twice its current
 size. Actually, just grow the array of aligned sequences.
 Copy the current aligned sequences into the new array
 free the now unused memory.
 Return 1 if success
 0 if failure
 */
int grow_alns_map_alignment(MapAlignmentP aln) {
    int i, j, k;
    int new_size;
    AlnSeqP as;
    AlnSeqP first_seq;
    AlnSeqP* NewAlnSeqArray;

    new_size = (aln->size) * 2;

    // Allocate another chunk of memory for AlnSeq[] as big as it
    // is now so we will double the size
    first_seq = (AlnSeqP) save_malloc(aln->size * sizeof (AlnSeq));
    if (first_seq == NULL) {
        fprintf(stderr, "Out of memory, sucka!\n");
        return 0;
    }

    // Now, allocate the new array of pointers to the aligned seqs
    NewAlnSeqArray = (AlnSeqP*) save_malloc(new_size * sizeof (AlnSeqP));
    if (NewAlnSeqArray == NULL) {
        fprintf(stderr, "Out of memory, sucka!\n");
        return 0;
    }

    // Now, point the pointers to the pointees
    // First, the old pointers/pointees
    for (i = 0; i < aln->size; i++) {
        NewAlnSeqArray[i] = aln->AlnSeqArray[i];
    }
    // Now, the new pointers/pointees
    k = 0;
    for (i = aln->size; i < new_size; i++) {
        /* Just in case there's some cruffy leftovers in our
         clean new memories */
        NewAlnSeqArray[i] = &first_seq[k++];
        /* Zero them out */
        as = NewAlnSeqArray[i];
        for (j = 0; j <= MAX_ID_LEN; j++) {
            as->id[j] = '\0';
        }
        for (j = 0; j <= MAX_DESC_LEN; j++) {
            as->desc[j] = '\0';
        }
        for (j = 0; j <= INIT_ALN_SEQ_LEN; j++) {
            as->seq[j] = '\0';
        }
        /* Set all their char* ins to NULL */
        for (j = 0; j <= INIT_ALN_SEQ_LEN; j++) {
            as->ins[j] = NULL;
        }
        as->start = 0;
        as->end = 0;
        as->revcom = 0;
        as->trimmed = 0;
        as->score = 0;
        as->segment = 'n';
    }

    // Now, the old aln->AlnSeqArray can be freed like a bird
    free(aln->AlnSeqArray);

    // And put in it's place the newer, bigger NewAlnSeqArray
    aln->AlnSeqArray = NewAlnSeqArray;
    aln->size = new_size;
    return 1;
}
