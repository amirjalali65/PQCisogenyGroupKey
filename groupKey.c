/******************************************************************************************************
 *      Supersingular Isogeny Grroup-Key implementation
 * 
 *  This file contains the main operations for performing an instance of supersingular isogeny 
 *  group key operations between three parties. The computation of 3- and 4- isogenies are 
 *  implemented based on the SIKE and SIDH implementation developed by Microsoft Research.
 *  This work is an extension to the SIDH library developed by Microsoft Research.
 *
 *  Modified and created by Amir Jalali             ajalali2016@fau.edu 
 ******************************************************************************************************/

#include "P747_internal.h"
#include "random.h"

static void clear_words(void *mem, digit_t nwords)
{   // Clear digits from memory. "nwords" indicates the number of digits to be zeroed.
    // This function uses the volatile type qualifier to inform the compiler not to optimize out the memory clearing.
    unsigned int i;
    volatile digit_t *v = mem;

    for (i = 0; i < nwords; i++)
    {
        v[i] = 0;
    }
}

static void init_basis(digit_t *gen, f2elm_t XP, f2elm_t XQ, f2elm_t XR)
{ // Initialization of basis points

    fpcopy(gen, XP[0]);
    fpcopy(gen + NWORDS_FIELD, XP[1]);
    fpcopy(gen + 2 * NWORDS_FIELD, XQ[0]);
    fpzero(XQ[1]);
    fpcopy(gen + 3 * NWORDS_FIELD, XR[0]);
    fpcopy(gen + 4 * NWORDS_FIELD, XR[1]);
}

static void init_alpha(digit_t *alpha, f2elm_t Alpha)
{ // Initialization of alpha aka 2-torsion point on the base curve
    fpzero(Alpha[0]);
    fpcopy(alpha, Alpha[1]);
}

static void fp2_encode(const f2elm_t x, unsigned char *enc)
{ // Conversion of GF(p^2) element from Montgomery to standard representation, and encoding by removing leading 0 bytes
    unsigned int i;
    f2elm_t t;

    from_fp2mont(x, t);
    for (i = 0; i < FP2_ENCODED_BYTES / 2; i++)
    {
        enc[i] = ((unsigned char *)t)[i];
        enc[i + FP2_ENCODED_BYTES / 2] = ((unsigned char *)t)[i + MAXBITS_FIELD / 8];
    }
}

static void fp2_decode(const unsigned char *enc, f2elm_t x)
{ // Parse byte sequence back into GF(p^2) element, and conversion to Montgomery representation
    unsigned int i;

    for (i = 0; i < 2 * (MAXBITS_FIELD / 8); i++)
        ((unsigned char *)x)[i] = 0;
    for (i = 0; i < FP2_ENCODED_BYTES / 2; i++)
    {
        ((unsigned char *)x)[i] = enc[i];
        ((unsigned char *)x)[i + MAXBITS_FIELD / 8] = enc[i + FP2_ENCODED_BYTES / 2];
    }
    to_fp2mont(x, x);
}

void random_mod_order_A(unsigned char *random_digits)
{   // Generation of Alice's secret key
    // Outputs random value in [0, 2^eA - 1]
    unsigned long long nbytes = NBITS_TO_NBYTES(OALICE_BITS);

    clear_words((void *)random_digits, MAXWORDS_ORDER);
    randombytes(random_digits, nbytes);
    random_digits[nbytes - 1] &= MASK_ALICE; // Masking last byte
}

void random_mod_order_B(unsigned char *random_digits)
{   // Generation of Bob's secret key
    // Outputs random value in [0, 2^Floor(Log(2, oB)) - 1]
    unsigned long long nbytes = NBITS_TO_NBYTES(OBOB_BITS - 1);

    clear_words((void *)random_digits, MAXWORDS_ORDER);
    randombytes(random_digits, nbytes);
    random_digits[nbytes - 1] &= MASK_BOB; // Masking last byte
}

void random_mod_order_C(unsigned char *random_digits)
{ // Generation of Eve's secret key
    // Outputs random value in [0, 2^Floor(Log(2, oC)) - 1]
    unsigned long long nbytes = NBITS_TO_NBYTES(OEVE_BITS - 1);

    clear_words((void *)random_digits, MAXWORDS_ORDER);
    randombytes(random_digits, nbytes);
    random_digits[nbytes - 1] &= MASK_EVE; // Masking last byte
}

int EphemeralKeyGeneration_A(const unsigned char *PrivateKeyA, unsigned char *PublicKeyA)
{   // Alice's ephemeral public key generation
    // Input:  a private key PrivateKeyA in the range [0, 2^eA - 1].
    // Output: the public key PublicKeyA consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, phiPB = {0}, phiQB = {0}, phiRB = {0}, phiPC = {0}, phiQC = {0}, phiRC = {0}, pts[MAX_INT_POINTS_ALICE];
    f2elm_t XPA, XQA, XRA, coeff[3], A24plus = {0}, C24 = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;

    // Initialize basis points
    init_basis((digit_t *)A_gen, XPA, XQA, XRA);
    init_basis((digit_t *)B_gen, phiPB->X, phiQB->X, phiRB->X);
    init_basis((digit_t *)C_gen, phiPC->X, phiQC->X, phiRC->X);
    //init_alpha((digit_t*)E0_alpha, phiAlpha->X);

    to_fp2mont(XPA, XPA);
    to_fp2mont(XQA, XQA);
    to_fp2mont(XRA, XRA);
    to_fp2mont(phiPB->X, phiPB->X);
    to_fp2mont(phiQB->X, phiQB->X);
    to_fp2mont(phiRB->X, phiRB->X);
    to_fp2mont(phiPC->X, phiPC->X);
    to_fp2mont(phiQC->X, phiQC->X);
    to_fp2mont(phiRC->X, phiRC->X);

    fpcopy((digit_t *)&Montgomery_one, (phiPB->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiQB->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiRB->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiPC->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiQC->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiRC->Z)[0]);

    // Initialize constants
    fpcopy((digit_t *)&Montgomery_one, A24plus[0]);
    fp2add(A24plus, A24plus, C24);

    // Retrieve kernel point
    LADDER3PT(XPA, XQA, XRA, (digit_t *)PrivateKeyA, ALICE, R, A);
    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Alice; row++)
    {
        while (index < MAX_Alice - row)
        {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Alice[MAX_Alice - index - row];
            xDBLe(R, R, A24plus, C24, (int)(2 * m));
            index += m;
        }

        get_4_isog(R, A24plus, C24, coeff);

        for (i = 0; i < npts; i++)
        {
            eval_4_isog(pts[i], coeff);
        }
        eval_4_isog(phiPB, coeff);
        eval_4_isog(phiQB, coeff);
        eval_4_isog(phiRB, coeff);
        eval_4_isog(phiPC, coeff);
        eval_4_isog(phiQC, coeff);
        eval_4_isog(phiRC, coeff);

        fp2copy(pts[npts - 1]->X, R->X);
        fp2copy(pts[npts - 1]->Z, R->Z);
        index = pts_index[npts - 1];
        npts -= 1;
    }

    get_4_isog(R, A24plus, C24, coeff);
    eval_4_isog(phiPB, coeff);
    eval_4_isog(phiQB, coeff);
    eval_4_isog(phiRB, coeff);
    eval_4_isog(phiPC, coeff);
    eval_4_isog(phiQC, coeff);
    eval_4_isog(phiRC, coeff);

    inv_6_way(phiPB->Z, phiQB->Z, phiRB->Z, phiPC->Z, phiQC->Z, phiRC->Z);
    fp2mul_mont(phiPB->X, phiPB->Z, phiPB->X);
    fp2mul_mont(phiQB->X, phiQB->Z, phiQB->X);
    fp2mul_mont(phiRB->X, phiRB->Z, phiRB->X);
    fp2mul_mont(phiPC->X, phiPC->Z, phiPC->X);
    fp2mul_mont(phiQC->X, phiQC->Z, phiQC->X);
    fp2mul_mont(phiRC->X, phiRC->Z, phiRC->X);

    // Format public key
    fp2_encode(phiPB->X, PublicKeyA);
    fp2_encode(phiQB->X, PublicKeyA + FP2_ENCODED_BYTES);
    fp2_encode(phiRB->X, PublicKeyA + 2 * FP2_ENCODED_BYTES);
    fp2_encode(phiPC->X, PublicKeyA + 3 * FP2_ENCODED_BYTES);
    fp2_encode(phiQC->X, PublicKeyA + 4 * FP2_ENCODED_BYTES);
    fp2_encode(phiRC->X, PublicKeyA + 5 * FP2_ENCODED_BYTES);

    return 0;
}

int EphemeralKeyGeneration_B(const unsigned char *PrivateKeyB, unsigned char *PublicKeyB)
{   // Bob's ephemeral public key generation
    // Input:  a private key PrivateKeyB in the range [0, 2^Floor(Log(2,oB)) - 1].
    // Output: the public key PublicKeyB consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, phiPA = {0}, phiQA = {0}, phiRA = {0}, phiPC = {0}, phiQC = {0}, phiRC = {0}, pts[MAX_INT_POINTS_BOB];
    f2elm_t XPB, XQB, XRB, coeff[3], A24plus = {0}, A24minus = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;

    // Initialize basis points
    init_basis((digit_t *)B_gen, XPB, XQB, XRB);
    init_basis((digit_t *)A_gen, phiPA->X, phiQA->X, phiRA->X);
    init_basis((digit_t *)C_gen, phiPC->X, phiQC->X, phiRC->X);

    to_fp2mont(XPB, XPB);
    to_fp2mont(XQB, XQB);
    to_fp2mont(XRB, XRB);
    to_fp2mont(phiPA->X, phiPA->X);
    to_fp2mont(phiQA->X, phiQA->X);
    to_fp2mont(phiRA->X, phiRA->X);
    to_fp2mont(phiPC->X, phiPC->X);
    to_fp2mont(phiQC->X, phiQC->X);
    to_fp2mont(phiRC->X, phiRC->X);

    fpcopy((digit_t *)&Montgomery_one, (phiPA->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiQA->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiRA->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiPC->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiQC->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiRC->Z)[0]);

    // Initialize constants
    fpcopy((digit_t *)&Montgomery_one, A24plus[0]);
    fp2add(A24plus, A24plus, A24plus);
    fp2copy(A24plus, A24minus);
    fp2neg(A24minus);

    // Retrieve kernel point
    LADDER3PT(XPB, XQB, XRB, (digit_t *)PrivateKeyB, BOB, R, A);
    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Bob; row++)
    {
        while (index < MAX_Bob - row)
        {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Bob[MAX_Bob - index - row];
            xTPLe(R, R, A24minus, A24plus, (int)m);
            index += m;
        }
        get_3_isog(R, A24minus, A24plus, coeff);

        for (i = 0; i < npts; i++)
        {
            eval_3_isog(pts[i], coeff);
        }
        eval_3_isog(phiPA, coeff);
        eval_3_isog(phiQA, coeff);
        eval_3_isog(phiRA, coeff);
        eval_3_isog(phiPC, coeff);
        eval_3_isog(phiQC, coeff);
        eval_3_isog(phiRC, coeff);

        fp2copy(pts[npts - 1]->X, R->X);
        fp2copy(pts[npts - 1]->Z, R->Z);
        index = pts_index[npts - 1];
        npts -= 1;
    }
    get_3_isog(R, A24minus, A24plus, coeff);
    eval_3_isog(phiPA, coeff);
    eval_3_isog(phiQA, coeff);
    eval_3_isog(phiRA, coeff);
    eval_3_isog(phiPC, coeff);
    eval_3_isog(phiQC, coeff);
    eval_3_isog(phiRC, coeff);

    inv_6_way(phiPA->Z, phiQA->Z, phiRA->Z, phiPC->Z, phiQC->Z, phiRC->Z);
    fp2mul_mont(phiPA->X, phiPA->Z, phiPA->X);
    fp2mul_mont(phiQA->X, phiQA->Z, phiQA->X);
    fp2mul_mont(phiRA->X, phiRA->Z, phiRA->X);
    fp2mul_mont(phiPC->X, phiPC->Z, phiPC->X);
    fp2mul_mont(phiQC->X, phiQC->Z, phiQC->X);
    fp2mul_mont(phiRC->X, phiRC->Z, phiRC->X);

    // Format public key
    fp2_encode(phiPA->X, PublicKeyB);
    fp2_encode(phiQA->X, PublicKeyB + FP2_ENCODED_BYTES);
    fp2_encode(phiRA->X, PublicKeyB + 2 * FP2_ENCODED_BYTES);
    fp2_encode(phiPC->X, PublicKeyB + 3 * FP2_ENCODED_BYTES);
    fp2_encode(phiQC->X, PublicKeyB + 4 * FP2_ENCODED_BYTES);
    fp2_encode(phiRC->X, PublicKeyB + 5 * FP2_ENCODED_BYTES);

    return 0;
}

int EphemeralKeyGeneration_C(const unsigned char *PrivateKeyC, unsigned char *PublicKeyC)
{   // Eve's ephemral public key generation
    // Input: a private key PrivateKeyC in the range [0, 2^Floor(Log(2,oC)) - 1].
    // Output: the public key PublicKeyC consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.

    // alpha is a point of order 2 on the base curve
    point_proj_t phiAlpha = {0}, R, R_2, phiPA = {0}, phiQA = {0}, phiRA = {0}, phiPB = {0}, phiQB = {0}, phiRB = {0}, pts[MAX_INT_POINTS_EVE];
    f2elm_t XPC, XQC, XRC, A24plus = {0}, C24 = {0}, A24pluscpy = {0}, C24cpy = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_EVE], npts = 0, ii = 0;

    // Initialized basis points
    init_basis((digit_t *)C_gen, XPC, XQC, XRC);
    init_basis((digit_t *)A_gen, phiPA->X, phiQA->X, phiRA->X);
    init_basis((digit_t *)B_gen, phiPB->X, phiQB->X, phiRB->X);
    init_alpha((digit_t *)E0_alpha, phiAlpha->X);

    to_fp2mont(phiAlpha->X, phiAlpha->X);
    to_fp2mont(XPC, XPC);
    to_fp2mont(XQC, XQC);
    to_fp2mont(XRC, XRC);
    to_fp2mont(phiPB->X, phiPB->X);
    to_fp2mont(phiQB->X, phiQB->X);
    to_fp2mont(phiRB->X, phiRB->X);
    to_fp2mont(phiPA->X, phiPA->X);
    to_fp2mont(phiQA->X, phiQA->X);
    to_fp2mont(phiRA->X, phiRA->X);

    fpcopy((digit_t *)&Montgomery_one, (phiPA->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiQA->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiRA->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiPB->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiQB->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiRB->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiAlpha->Z)[0]);

    // Initialize constants
    fpcopy((digit_t *)&Montgomery_one, A24plus[0]);
    fp2add(A24plus, A24plus, C24);

    // Retrieve kernel point
    LADDER3PT(XPC, XQC, XRC, (digit_t *)PrivateKeyC, EVE, R, A);

    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Eve; row++)
    {
        while (index < MAX_Eve - row)
        {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Eve[MAX_Eve - index - row];
            xQNTPLe(R, R, A24plus, C24, (int)m);
            index += m;
        }
        xDBL(R, R_2, A24plus, C24);
        eval_5_isog(R, R_2, phiAlpha);

        for (i = 0; i < npts; i++)
        {
            eval_5_isog(R, R_2, pts[i]);
        }
        eval_5_isog(R, R_2, phiPA);
        eval_5_isog(R, R_2, phiQA);
        eval_5_isog(R, R_2, phiRA);
        eval_5_isog(R, R_2, phiPB);
        eval_5_isog(R, R_2, phiQB);
        eval_5_isog(R, R_2, phiRB);
        get_a_from_alpha(phiAlpha, A24plus, C24);

        fp2copy(pts[npts - 1]->X, R->X);
        fp2copy(pts[npts - 1]->Z, R->Z);
        index = pts_index[npts - 1];
        npts -= 1;
    }
    xDBL(R, R_2, A24plus, C24);
    eval_5_isog(R, R_2, phiPA);
    eval_5_isog(R, R_2, phiQA);
    eval_5_isog(R, R_2, phiRA);
    eval_5_isog(R, R_2, phiPB);
    eval_5_isog(R, R_2, phiQB);
    eval_5_isog(R, R_2, phiRB);

    inv_6_way(phiPA->Z, phiQA->Z, phiRA->Z, phiPB->Z, phiQB->Z, phiRB->Z);
    fp2mul_mont(phiPA->X, phiPA->Z, phiPA->X);
    fp2mul_mont(phiQA->X, phiQA->Z, phiQA->X);
    fp2mul_mont(phiRA->X, phiRA->Z, phiRA->X);
    fp2mul_mont(phiPB->X, phiPB->Z, phiPB->X);
    fp2mul_mont(phiQB->X, phiQB->Z, phiQB->X);
    fp2mul_mont(phiRB->X, phiRB->Z, phiRB->X);

    // Format public key
    fp2_encode(phiPA->X, PublicKeyC);
    fp2_encode(phiQA->X, PublicKeyC + FP2_ENCODED_BYTES);
    fp2_encode(phiRA->X, PublicKeyC + 2 * FP2_ENCODED_BYTES);
    fp2_encode(phiPB->X, PublicKeyC + 3 * FP2_ENCODED_BYTES);
    fp2_encode(phiQB->X, PublicKeyC + 4 * FP2_ENCODED_BYTES);
    fp2_encode(phiRB->X, PublicKeyC + 5 * FP2_ENCODED_BYTES);

    return 0;
}

int BSharedPublicFromA(const unsigned char *PrivateKeyB, const unsigned char *PublicKeyA, unsigned char *SharedPublicAB)
{   // Bob's ephemeral shared secret computation
    // It produces a shared secret key SharedSecretB using his secret key PrivateKeyB and Alice's public key PublicKeyA
    // Inputs: Bob's PrivateKeyB is an integer in the range [0, 2^Floor(Log(2,oB)) - 1].
    //         Alice's PublicKeyA consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
    // Output: a shared secret SharedSecretB that consists of one element in GF(p^2) encoded by removing leading 0 bytes.
    point_proj_t R, phiA_PC = {0}, phiA_QC = {0}, phiA_RC = {0}, pts[MAX_INT_POINTS_BOB];
    f2elm_t coeff[3], PKB[3];
    f2elm_t A24plus = {0}, A24minus = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;

    // Initialize images of Alice's basis
    fp2_decode(PublicKeyA, PKB[0]);
    fp2_decode(PublicKeyA + FP2_ENCODED_BYTES, PKB[1]);
    fp2_decode(PublicKeyA + 2 * FP2_ENCODED_BYTES, PKB[2]);
    fp2_decode(PublicKeyA + 3 * FP2_ENCODED_BYTES, phiA_PC->X);
    fp2_decode(PublicKeyA + 4 * FP2_ENCODED_BYTES, phiA_QC->X);
    fp2_decode(PublicKeyA + 5 * FP2_ENCODED_BYTES, phiA_RC->X);
    fpcopy((digit_t *)&Montgomery_one, (phiA_PC->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiA_QC->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiA_RC->Z)[0]);

    // Initialize constants
    get_A(PKB[0], PKB[1], PKB[2], A); // Retrieve E_A
    fpadd((digit_t *)&Montgomery_one, (digit_t *)&Montgomery_one, A24minus[0]);
    fp2add(A, A24minus, A24plus);
    fp2sub(A, A24minus, A24minus);

    // Retrieve kernel point
    LADDER3PT(PKB[0], PKB[1], PKB[2], (digit_t *)PrivateKeyB, BOB, R, A);

    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Bob; row++)
    {
        while (index < MAX_Bob - row)
        {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Bob[MAX_Bob - index - row];
            xTPLe(R, R, A24minus, A24plus, (int)m);
            index += m;
        }
        get_3_isog(R, A24minus, A24plus, coeff);

        for (i = 0; i < npts; i++)
        {
            eval_3_isog(pts[i], coeff);
        }

        eval_3_isog(phiA_PC, coeff);
        eval_3_isog(phiA_QC, coeff);
        eval_3_isog(phiA_RC, coeff);

        fp2copy(pts[npts - 1]->X, R->X);
        fp2copy(pts[npts - 1]->Z, R->Z);
        index = pts_index[npts - 1];
        npts -= 1;
    }
    get_3_isog(R, A24minus, A24plus, coeff);
    eval_3_isog(phiA_PC, coeff);
    eval_3_isog(phiA_QC, coeff);
    eval_3_isog(phiA_RC, coeff);

    inv_3_way(phiA_PC->Z, phiA_QC->Z, phiA_RC->Z);
    fp2mul_mont(phiA_PC->X, phiA_PC->Z, phiA_PC->X);
    fp2mul_mont(phiA_QC->X, phiA_QC->Z, phiA_QC->X);
    fp2mul_mont(phiA_RC->X, phiA_RC->Z, phiA_RC->X);

    fp2_encode(phiA_PC->X, SharedPublicAB);                         // phi_AB(PC)
    fp2_encode(phiA_QC->X, SharedPublicAB + FP2_ENCODED_BYTES);     // phi_AB(QC)
    fp2_encode(phiA_RC->X, SharedPublicAB + 2 * FP2_ENCODED_BYTES); // phi_AB(RC)

    return 0;
}

int CSharedSecretFromB(const unsigned char *PrivateKeyC, const unsigned char *PublicKeyB, const unsigned char *SharedPublicAB, unsigned char *SharedPublicBC, unsigned char *SharedSecret)
{
    point_proj_t R, R_2 = {0}, phiAB_PC = {0}, phiAB_QC = {0}, phiAB_RC = {0}, phiB_PA = {0}, phiB_QA = {0}, phiB_RA = {0}, pts[MAX_INT_POINTS_EVE];
    f2elm_t PKB[3], PKAB[3], jinv;
    f2elm_t A24plus = {0}, C24 = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_EVE], npts = 0, ii = 0;

    // Initialize images from B
    fp2_decode(PublicKeyB, phiB_PA->X);
    fp2_decode(PublicKeyB + FP2_ENCODED_BYTES, phiB_QA->X);
    fp2_decode(PublicKeyB + 2 * FP2_ENCODED_BYTES, phiB_RA->X);
    fp2_decode(PublicKeyB + 3 * FP2_ENCODED_BYTES, PKB[0]);
    fp2_decode(PublicKeyB + 4 * FP2_ENCODED_BYTES, PKB[1]);
    fp2_decode(PublicKeyB + 5 * FP2_ENCODED_BYTES, PKB[2]);

    fpcopy((digit_t *)&Montgomery_one, (phiB_PA->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiB_QA->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiB_RA->Z)[0]);

    // Initialize images from AB
    fp2_decode(SharedPublicAB, PKAB[0]);
    fp2_decode(SharedPublicAB + FP2_ENCODED_BYTES, PKAB[1]);
    fp2_decode(SharedPublicAB + 2 * FP2_ENCODED_BYTES, PKAB[2]);

    get_A(PKB[0], PKB[1], PKB[2], A); // Retrieve E_B
    fpadd((digit_t *)&Montgomery_one, (digit_t *)&Montgomery_one, C24[0]);
    fp2add(A, C24, A24plus);
    fpadd(C24[0], C24[0], C24[0]);

    // Retrieve kernel point
    LADDER3PT(PKB[0], PKB[1], PKB[2], (digit_t *)PrivateKeyC, EVE, R, A);

    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Eve; row++)
    {
        while (index < MAX_Eve - row)
        {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Eve[MAX_Eve - index - row];
            xQNTPLe(R, R, A24plus, C24, (int)m);
            index += m;
        }
        xDBL(R, R_2, A24plus, C24);
        eval_5_isog(R, R_2, phiB_PA);
        eval_5_isog(R, R_2, phiB_QA);
        eval_5_isog(R, R_2, phiB_RA);

        get_A_projective(phiB_PA, phiB_QA, phiB_RA, A24plus, C24);

        for (i = 0; i < npts; i++)
        {
            eval_5_isog(R, R_2, pts[i]);
        }

        fp2copy(pts[npts - 1]->X, R->X);
        fp2copy(pts[npts - 1]->Z, R->Z);
        index = pts_index[npts - 1];
        npts -= 1;
    }
    xDBL(R, R_2, A24plus, C24);
    eval_5_isog(R, R_2, phiB_PA);
    eval_5_isog(R, R_2, phiB_QA);
    eval_5_isog(R, R_2, phiB_RA);

    inv_3_way(phiB_PA->Z, phiB_QA->Z, phiB_RA->Z);
    fp2mul_mont(phiB_PA->X, phiB_PA->Z, phiB_PA->X);
    fp2mul_mont(phiB_QA->X, phiB_QA->Z, phiB_QA->X);
    fp2mul_mont(phiB_RA->X, phiB_RA->Z, phiB_RA->X);

    // Format public key
    fp2_encode(phiB_PA->X, SharedPublicBC);                         // phic(phiB(PA))
    fp2_encode(phiB_QA->X, SharedPublicBC + FP2_ENCODED_BYTES);     // phic(phiB(QA))
    fp2_encode(phiB_RA->X, SharedPublicBC + 2 * FP2_ENCODED_BYTES); // phic(phiB(RA))

    fp2zero(A);
    fp2zero(A24plus);
    fp2zero(C24);
    // Computing the shared secret
    get_A(PKAB[0], PKAB[1], PKAB[2], A); // Retrieve E_AB
    fpadd((digit_t *)&Montgomery_one, (digit_t *)&Montgomery_one, C24[0]);
    fp2add(A, C24, A24plus);
    fpadd(C24[0], C24[0], C24[0]);

    // Retrieve kernel point
    LADDER3PT(PKAB[0], PKAB[1], PKAB[2], (digit_t *)PrivateKeyC, EVE, R, A);

    // Initialize the phiAB points
    fp2copy(PKAB[0], phiAB_PC->X);
    fp2copy(PKAB[1], phiAB_QC->X);
    fp2copy(PKAB[2], phiAB_RC->X);

    fpcopy((digit_t *)&Montgomery_one, (phiAB_PC->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiAB_QC->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiAB_RC->Z)[0]);

    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Eve; row++)
    {
        while (index < MAX_Eve - row)
        {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Eve[MAX_Eve - index - row];
            xQNTPLe(R, R, A24plus, C24, (int)(m));
            index += m;
        }
        xDBL(R, R_2, A24plus, C24);
        eval_5_isog(R, R_2, phiAB_PC);
        eval_5_isog(R, R_2, phiAB_QC);
        eval_5_isog(R, R_2, phiAB_RC);
        get_A_projective(phiAB_PC, phiAB_QC, phiAB_RC, A24plus, C24);

        for (i = 0; i < npts; i++)
        {
            eval_5_isog(R, R_2, pts[i]);
        }

        fp2copy(pts[npts - 1]->X, R->X);
        fp2copy(pts[npts - 1]->Z, R->Z);
        index = pts_index[npts - 1];
        npts -= 1;
    }
    xDBL(R, R_2, A24plus, C24);
    eval_5_isog(R, R_2, phiAB_PC);
    eval_5_isog(R, R_2, phiAB_QC);
    eval_5_isog(R, R_2, phiAB_RC);
    get_A_projective(phiAB_PC, phiAB_QC, phiAB_RC, A24plus, C24);

    fp2div2(C24, C24);
    fp2sub(A24plus, C24, A24plus);
    fp2div2(C24, C24);
    j_inv(A24plus, C24, jinv);
    fp2_encode(jinv, SharedSecret); // Format shared secret

    return 0;
}

int ASharedSecretFromC(const unsigned char *PrivateKeyA, const unsigned char *PublicKeyC, const unsigned char *SharedPublicBC, unsigned char *SharedPublicAC, unsigned char *SharedSecret)
{
    point_proj_t R, phiC_PB = {0}, phiC_QB = {0}, phiC_RB = {0}, pts[MAX_INT_POINTS_EVE];
    f2elm_t coeff[3], PKC[3], PKBC[3], jinv;
    f2elm_t A24plus = {0}, C24 = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_EVE], npts = 0, ii = 0;

    // Initialize images from C
    fp2_decode(PublicKeyC, PKC[0]);
    fp2_decode(PublicKeyC + FP2_ENCODED_BYTES, PKC[1]);
    fp2_decode(PublicKeyC + 2 * FP2_ENCODED_BYTES, PKC[2]);
    fp2_decode(PublicKeyC + 3 * FP2_ENCODED_BYTES, phiC_PB->X);
    fp2_decode(PublicKeyC + 4 * FP2_ENCODED_BYTES, phiC_QB->X);
    fp2_decode(PublicKeyC + 5 * FP2_ENCODED_BYTES, phiC_RB->X);

    fpcopy((digit_t *)&Montgomery_one, (phiC_PB->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiC_QB->Z)[0]);
    fpcopy((digit_t *)&Montgomery_one, (phiC_RB->Z)[0]);

    // Initialize images from E_AB
    fp2_decode(SharedPublicBC, PKBC[0]);
    fp2_decode(SharedPublicBC + FP2_ENCODED_BYTES, PKBC[1]);
    fp2_decode(SharedPublicBC + 2 * FP2_ENCODED_BYTES, PKBC[2]);

    get_A(PKC[0], PKC[1], PKC[2], A); // Retrieve E_C
    fpadd((digit_t *)&Montgomery_one, (digit_t *)&Montgomery_one, C24[0]);
    fp2add(A, C24, A24plus);
    fpadd(C24[0], C24[0], C24[0]);

    // Retrieve kernel point
    LADDER3PT(PKC[0], PKC[1], PKC[2], (digit_t *)PrivateKeyA, ALICE, R, A);

    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Alice; row++)
    {
        while (index < MAX_Alice - row)
        {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Alice[MAX_Alice - index - row];
            xDBLe(R, R, A24plus, C24, (int)(2 * m));
            index += m;
        }
        get_4_isog(R, A24plus, C24, coeff);

        for (i = 0; i < npts; i++)
        {
            eval_4_isog(pts[i], coeff);
        }
        eval_4_isog(phiC_PB, coeff);
        eval_4_isog(phiC_QB, coeff);
        eval_4_isog(phiC_RB, coeff);

        fp2copy(pts[npts - 1]->X, R->X);
        fp2copy(pts[npts - 1]->Z, R->Z);
        index = pts_index[npts - 1];
        npts -= 1;
    }
    get_4_isog(R, A24plus, C24, coeff);
    eval_4_isog(phiC_PB, coeff);
    eval_4_isog(phiC_QB, coeff);
    eval_4_isog(phiC_RB, coeff);

    inv_3_way(phiC_PB->Z, phiC_QB->Z, phiC_RB->Z);
    fp2mul_mont(phiC_PB->X, phiC_PB->Z, phiC_PB->X);
    fp2mul_mont(phiC_QB->X, phiC_QB->Z, phiC_QB->X);
    fp2mul_mont(phiC_RB->X, phiC_RB->Z, phiC_RB->X);

    // Format public key
    fp2_encode(phiC_PB->X, SharedPublicAC);                         // phiA(phiC(PB))
    fp2_encode(phiC_QB->X, SharedPublicAC + FP2_ENCODED_BYTES);     // phiA(phiC(QB))
    fp2_encode(phiC_RB->X, SharedPublicAC + 2 * FP2_ENCODED_BYTES); // phiA(phiC(RB))

    fp2zero(A);
    fp2zero(A24plus);
    fp2zero(C24);
    // Compute the shared secret
    get_A(PKBC[0], PKBC[1], PKBC[2], A); // Retrieve E_C
    fpadd((digit_t *)&Montgomery_one, (digit_t *)&Montgomery_one, C24[0]);
    fp2add(A, C24, A24plus);
    fpadd(C24[0], C24[0], C24[0]);

    // Retrieve kernel point
    LADDER3PT(PKBC[0], PKBC[1], PKBC[2], (digit_t *)PrivateKeyA, ALICE, R, A);

    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Alice; row++)
    {
        while (index < MAX_Alice - row)
        {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Alice[MAX_Alice - index - row];
            xDBLe(R, R, A24plus, C24, (int)(2 * m));
            index += m;
        }
        get_4_isog(R, A24plus, C24, coeff);

        for (i = 0; i < npts; i++)
        {
            eval_4_isog(pts[i], coeff);
        }

        fp2copy(pts[npts - 1]->X, R->X);
        fp2copy(pts[npts - 1]->Z, R->Z);
        index = pts_index[npts - 1];
        npts -= 1;
    }
    get_4_isog(R, A24plus, C24, coeff);
    fp2div2(C24, C24);
    fp2sub(A24plus, C24, A24plus);
    fp2div2(C24, C24);
    j_inv(A24plus, C24, jinv);
    fp2_encode(jinv, SharedSecret); // Format shared secret

    return 0;
}

int BSharedSecretFromA(const unsigned char *PrivateKeyB, const unsigned char *SharedPublicAC, unsigned char *SharedSecret)
{
    point_proj_t R, pts[MAX_INT_POINTS_BOB];
    f2elm_t coeff[3], PKAC[3], jinv;
    f2elm_t A24plus = {0}, A24minus = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;

    // Initialize images of Alice's basis
    fp2_decode(SharedPublicAC, PKAC[0]);
    fp2_decode(SharedPublicAC + FP2_ENCODED_BYTES, PKAC[1]);
    fp2_decode(SharedPublicAC + 2 * FP2_ENCODED_BYTES, PKAC[2]);

    // Initialize constants
    get_A(PKAC[0], PKAC[1], PKAC[2], A); // Retrieve E_AC
    fpadd((digit_t *)&Montgomery_one, (digit_t *)&Montgomery_one, A24minus[0]);
    fp2add(A, A24minus, A24plus);
    fp2sub(A, A24minus, A24minus);

    // Retrieve kernel point
    LADDER3PT(PKAC[0], PKAC[1], PKAC[2], (digit_t *)PrivateKeyB, BOB, R, A);

    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Bob; row++)
    {
        while (index < MAX_Bob - row)
        {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Bob[MAX_Bob - index - row];
            xTPLe(R, R, A24minus, A24plus, (int)m);
            index += m;
        }
        get_3_isog(R, A24minus, A24plus, coeff);

        for (i = 0; i < npts; i++)
        {
            eval_3_isog(pts[i], coeff);
        }

        fp2copy(pts[npts - 1]->X, R->X);
        fp2copy(pts[npts - 1]->Z, R->Z);
        index = pts_index[npts - 1];
        npts -= 1;
    }
    get_3_isog(R, A24minus, A24plus, coeff);
    fp2add(A24plus, A24minus, A);
    fp2add(A, A, A);
    fp2sub(A24plus, A24minus, A24plus);
    j_inv(A, A24plus, jinv);
    fp2_encode(jinv, SharedSecret); // Format shared secret

    return 0;
}
