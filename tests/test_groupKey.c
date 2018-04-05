/********************************************************************************************
* Supersingular Isogeny Group Key Agreement Library
*
* Abstract: benchmarking/testing isogeny-based group key agreement mechanism
*********************************************************************************************/ 

// Benchmark and test parameters 
#define BENCH_LOOPS       5       
#define TEST_LOOPS        5 


int cryptotest_groupkey()
{ // Testing GROUPKEY
    unsigned int i,j;
	unsigned char sk_A[CRYPTO_SECRETKEYBYTES] = { 0 };
	unsigned char sk_B[CRYPTO_SECRETKEYBYTES] = { 0 };
	unsigned char sk_C[CRYPTO_SECRETKEYBYTES] = { 0 };
	unsigned char pk_A[CRYPTO_PUBLICKEYBYTES] = { 0 };
	unsigned char pk_B[CRYPTO_PUBLICKEYBYTES] = { 0 };
	unsigned char pk_C[CRYPTO_PUBLICKEYBYTES] = { 0 };
	unsigned char sp_AB[CRYPTO_SHAREDPUBLICBYTES] = { 0 };
	unsigned char sp_BC[CRYPTO_SHAREDPUBLICBYTES] = { 0 };
	unsigned char sp_AC[CRYPTO_SHAREDPUBLICBYTES] = { 0 };
	unsigned char ss_A[CRYPTO_BYTES] = { 0 };
	unsigned char ss_B[CRYPTO_BYTES] = { 0 };
	unsigned char ss_C[CRYPTO_BYTES] = { 0 };
    bool passed = true;

    printf("\n\nTESTING ISOGENY-BASED GROUP KEY MECHANISM %s\n", SCHEME_NAME);
    printf("--------------------------------------------------------------------------------------------------------\n\n");

	for (i = 0; i < TEST_LOOPS; i++) 
    {
        // secret-key generation
        random_mod_order_A(sk_A);
        random_mod_order_B(sk_B);
        random_mod_order_C(sk_C);

        // public-key generation
        EphemeralKeyGeneration_A(sk_A, pk_A);
        EphemeralKeyGeneration_B(sk_B, pk_B);
        EphemeralKeyGeneration_C(sk_C, pk_C);

        // key-agreement
        BSharedPublicFromA(sk_B, pk_A, sp_AB);
        CSharedSecretFromB(sk_C, pk_B, sp_AB, sp_BC, ss_C);
        ASharedSecretFromC(sk_A, pk_C, sp_BC, sp_AC, ss_A);
        BSharedSecretFromA(sk_B, sp_AC, ss_B);
 
		if ((memcmp(ss_A, ss_B, CRYPTO_BYTES) || memcmp(ss_A, ss_C, CRYPTO_BYTES)) != 0) {
            passed = false;
            //break;
        }
    }

    if (passed == true) printf("  GROUP KEY tests .................................................... PASSED");
    else { printf("  GROUP KEY tests ... FAILED"); printf("\n"); return FAILED; }
    printf("\n"); 

    return PASSED;
}

int cryptorun_groupkey()
{ // Benchmarking group key exchange
    unsigned int n;
    unsigned char sk_A[CRYPTO_SECRETKEYBYTES] = {0};
    unsigned char sk_B[CRYPTO_SECRETKEYBYTES] = {0};
    unsigned char sk_C[CRYPTO_SECRETKEYBYTES] = {0};
    unsigned char pk_A[CRYPTO_PUBLICKEYBYTES] = {0};
    unsigned char pk_B[CRYPTO_PUBLICKEYBYTES] = {0};
    unsigned char pk_C[CRYPTO_PUBLICKEYBYTES] = {0};
    unsigned char sp_AB[CRYPTO_SHAREDPUBLICBYTES] = {0};
    unsigned char sp_BC[CRYPTO_SHAREDPUBLICBYTES] = {0};
    unsigned char sp_AC[CRYPTO_SHAREDPUBLICBYTES] = {0};    
    unsigned char ss_A[CRYPTO_BYTES] = {0};
    unsigned char ss_B[CRYPTO_BYTES] = {0};
    unsigned char ss_C[CRYPTO_BYTES] = {0};
    unsigned long long cycles, cycles1, cycles2;

    printf("\n\nBENCHMARKING ISOGENY-BASED GROUP KEY MECHANISM %s\n", SCHEME_NAME);
    printf("--------------------------------------------------------------------------------------------------------\n\n");

    // Benchmarking key generation
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        EphemeralKeyGeneration_A(sk_A, pk_A);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Key generation A runs in ....................................... %10lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        EphemeralKeyGeneration_B(sk_B, pk_B);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Key generation B runs in ....................................... %10lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        EphemeralKeyGeneration_C(sk_C, pk_C);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Key generation C runs in ....................................... %10lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // Benchmarking key agreement
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        BSharedPublicFromA(sk_B, pk_A, sp_AB);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  B sharedPublic from A runs in .................................. %10lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        CSharedSecretFromB(sk_C, pk_B, sp_AB, sp_BC, ss_C);  
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  C sharedSecret from B runs in .................................. %10lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        ASharedSecretFromC(sk_A, pk_C, sp_BC, sp_AC, ss_A);  
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  A sharedSecret from C runs in .................................. %10lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        BSharedSecretFromA(sk_B, sp_AC, ss_B);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  B sharedSecret from A runs in .................................. %10lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    return PASSED;
}

int main()
{
    int Status = PASSED;
    
    Status = cryptotest_groupkey();             // Test group key agreement
    if (Status != PASSED) {
        printf("\n\n   Error detected: GROUPKEY_ERROR_SHARED_KEY \n\n");
        //return FAILED;
    }

    Status = cryptorun_groupkey();              // Benchmark group key agreement
    if (Status != PASSED) {
        printf("\n\n   Error detected: GROUPKEY_ERROR_SHARED_KEY \n\n");
        //return FAILED;
    }

    return Status;
}
