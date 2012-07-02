/* 
 * File:   mia_testsuite.c
 * Author: michael_siebauer
 *
 * Created on April 2, 2009, 1:13 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <CUnit/Basic.h>
#include "../src/mia.h"
#include "../src/map_align.h"
#include "../src/map_alignment.h"
#include "../src/io.h"
#include "../src/fsdb.h"

 RefSeqP ref_seq;
 FragSeqP frag_seq;
 FSDB frag_db;


    int maxi(int i1, int i2)
    {
      return (i1 > i2) ? i1 : i2;
    }

    void test_maxi(void)
    {
      CU_ASSERT(maxi(0,2) == 2);
      CU_ASSERT(maxi(0,-2) == 0);
      CU_ASSERT(maxi(2,2) == 2);
    }

    int init_testsuite(void){
        ref_seq = (RefSeqP)calloc(1, sizeof(RefSeqP));
        frag_seq = (FragSeqP)calloc(1, sizeof(FragSeqP));
        frag_db = init_FSDB();


        // read in our test reference sequence
        if (read_fasta_ref(ref_seq, "tr1.fna") != 1)
            return EXIT_FAILURE;

        FILE* frag_file = fileOpen("tf.fna", "r");
        if (frag_file == NULL)
            return EXIT_FAILURE;

        while (read_fasta(frag_file, frag_seq)){
            printf("%s\n", frag_seq->id);
        }

        

        return EXIT_SUCCESS;
    }

    int cleanup_testsuite(void){
        free(ref_seq);
        free(frag_db);
        return EXIT_SUCCESS;
    }

/*
 * 
 */
int main(int argc, char** argv) {
    // handle for testsuite
    CU_pSuite test = NULL;
    

    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    // create the suite
    test = CU_add_suite("Testsuite 1", init_testsuite, cleanup_testsuite);
    if (NULL == test) {
      CU_cleanup_registry();
      return CU_get_error();
   }

    // add tests to suite
    CU_add_test(test, "Basic Test", test_maxi);


    // Now Run all tests
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();
}

