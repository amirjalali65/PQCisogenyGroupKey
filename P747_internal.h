/********************************************************************************************
* Supersingular Isogeny Group Key Agreement Library
*
* Abstract: internal header file for P747
* This file is created by modifying the P751_internal.h file in SIKE package developed by Microsoft Research
*
*   Modified by Amir Jalali             ajalali2016@fau.edu
*********************************************************************************************/

#ifndef __P747_INTERNAL_H__
#define __P747_INTERNAL_H__

#include "api.h"

#if (TARGET == TARGET_AMD64)
#define NWORDS_FIELD 12   // Number of words of a 747-bit field element
#define p747_ZERO_WORDS 4 // Number of "0" digits in the least significant part of p747 + 1
#elif (TARGET == TARGET_x86)
#define NWORDS_FIELD 24
#define p747_ZERO_WORDS 8
#elif (TARGET == TARGET_ARM)
#define NWORDS_FIELD 24
#define p747_ZERO_WORDS 8
#elif (TARGET == TARGET_ARM64)
#define NWORDS_FIELD 12
#define p747_ZERO_WORDS 4
#endif

// Basic constants

#define NBITS_FIELD 747
#define MAXBITS_FIELD 768
#define MAXWORDS_FIELD ((MAXBITS_FIELD + RADIX - 1) / RADIX) // Max. number of words to represent field elements
#define NWORDS64_FIELD ((NBITS_FIELD + 63) / 64)             // Number of 64-bit words of a 747-bit field element
#define NBITS_ORDER 320
#define NWORDS_ORDER ((NBITS_ORDER + RADIX - 1) / RADIX) // Number of words of oA and oB, where oA and oB are the subgroup orders of Alice and Bob, resp.
#define NWORDS64_ORDER ((NBITS_ORDER + 63) / 64)         // Number of 64-bit words of a 384-bit element
#define MAXBITS_ORDER NBITS_ORDER
#define MAXWORDS_ORDER ((MAXBITS_ORDER + RADIX - 1) / RADIX) // Max. number of words to represent elements in [1, oA-1] or [1, oB].
#define ALICE 0
#define BOB 1
#define EVE 2
#define OALICE_BITS 261
#define OBOB_BITS 243
#define OEVE_BITS 244
#define OBOB_EXPON 153
#define OEVE_EXPON 105
#define MASK_ALICE 0x00
#define MASK_BOB 0x00
#define MASK_EVE 0x00
#define PRIME p747
#define PARAM_A 0
#define PARAM_C 1
// Fixed parameters for isogeny tree computation
#define MAX_INT_POINTS_ALICE 8
#define MAX_INT_POINTS_BOB 10
#define MAX_INT_POINTS_EVE 11
#define MAX_Alice 130
#define MAX_Bob 153
#define MAX_Eve 105
#define SECRETKEY_A_BYTES (OALICE_BITS + 7) / 8
#define SECRETKEY_B_BYTES (OBOB_BITS + 7) / 8
#define SECRETKEY_E_BYTES (OEVE_BITS + 7) / 8
#define FP2_ENCODED_BYTES 2 * ((NBITS_FIELD + 7) / 8)

// SIDH's basic element definitions and point representations

typedef digit_t felm_t[NWORDS_FIELD];      // Datatype for representing 747-bit field elements (768-bit max.)
typedef digit_t dfelm_t[2 * NWORDS_FIELD]; // Datatype for representing double-precision 2x747-bit field elements (2x768-bit max.)
typedef felm_t f2elm_t[2];                 // Datatype for representing quadratic extension field elements GF(p747^2)
typedef f2elm_t publickey_t[3];            // Datatype for representing public keys equivalent to three GF(p747^2) elements

typedef struct
{
    f2elm_t X;
    f2elm_t Z;
} point_proj; // Point representation in projective XZ Montgomery coordinates.
typedef point_proj point_proj_t[1];

/**************** Function prototypes ****************/
/************* Multiprecision functions **************/

// Copy wordsize digits, c = a, where lng(a) = nwords
void copy_words(const digit_t *a, digit_t *c, const unsigned int nwords);

// Multiprecision addition, c = a+b, where lng(a) = lng(b) = nwords. Returns the carry bit
unsigned int mp_add(const digit_t *a, const digit_t *b, digit_t *c, const unsigned int nwords);

// 747-bit multiprecision addition, c = a+b
void mp_add747(const digit_t *a, const digit_t *b, digit_t *c);
void mp_add747_asm(const digit_t *a, const digit_t *b, digit_t *c);
//void mp_addmask747_asm(const digit_t* a, const digit_t mask, digit_t* c);

// 2x747-bit multiprecision addition, c = a+b
void mp_add747x2(const digit_t *a, const digit_t *b, digit_t *c);
void mp_add747x2_asm(const digit_t *a, const digit_t *b, digit_t *c);

// Multiprecision subtraction, c = a-b, where lng(a) = lng(b) = nwords. Returns the borrow bit
unsigned int mp_sub(const digit_t *a, const digit_t *b, digit_t *c, const unsigned int nwords);
digit_t mp_sub747x2_asm(const digit_t *a, const digit_t *b, digit_t *c);

// Multiprecision left shift
void mp_shiftleft(digit_t *x, unsigned int shift, const unsigned int nwords);

// Multiprecision right shift by one
void mp_shiftr1(digit_t *x, const unsigned int nwords);

// Multiprecision left right shift by one
void mp_shiftl1(digit_t *x, const unsigned int nwords);

// Digit multiplication, digit * digit -> 2-digit result
void digit_x_digit(const digit_t a, const digit_t b, digit_t *c);
void MUL64(const digit_t a, digit_t b, digit_t* c);

// Multiprecision comba multiply, c = a*b, where lng(a) = lng(b) = nwords.
void mp_mul(const digit_t *a, const digit_t *b, digit_t *c, const unsigned int nwords);

void multiply(const digit_t *a, const digit_t *b, digit_t *c, const unsigned int nwords);

// Montgomery multiplication modulo the group order, mc = ma*mb*r' mod order, where ma,mb,mc in [0, order-1]
void Montgomery_multiply_mod_order(const digit_t *ma, const digit_t *mb, digit_t *mc, const digit_t *order, const digit_t *Montgomery_rprime);

// (Non-constant time) Montgomery inversion modulo the curve order using a^(-1) = a^(order-2) mod order
//void Montgomery_inversion_mod_order(const digit_t* ma, digit_t* mc, const digit_t* order, const digit_t* Montgomery_rprime);

void Montgomery_inversion_mod_order_bingcd(const digit_t *a, digit_t *c, const digit_t *order, const digit_t *Montgomery_rprime, const digit_t *Montgomery_R2);

// Conversion of elements in Z_r to Montgomery representation, where the order r is up to 384 bits.
void to_Montgomery_mod_order(const digit_t *a, digit_t *mc, const digit_t *order, const digit_t *Montgomery_rprime, const digit_t *Montgomery_Rprime);

// Conversion of elements in Z_r from Montgomery to standard representation, where the order is up to 384 bits.
void from_Montgomery_mod_order(const digit_t *ma, digit_t *c, const digit_t *order, const digit_t *Montgomery_rprime);

// Inversion modulo Alice's order 2^372.
void inv_mod_orderA(const digit_t *a, digit_t *c);

/************ Field arithmetic functions *************/

// Copy of a field element, c = a
void fpcopy747(const felm_t a, felm_t c);

// Zeroing a field element, a = 0
void fpzero747(felm_t a);

// Non constant-time comparison of two field elements. If a = b return TRUE, otherwise, return FALSE
bool fpequal747_non_constant_time(const felm_t a, const felm_t b);

// Modular addition, c = a+b mod p747
extern void fpadd747(const digit_t *a, const digit_t *b, digit_t *c);
extern void fpadd747_asm(const digit_t *a, const digit_t *b, digit_t *c);

// Modular subtraction, c = a-b mod p747
extern void fpsub747(const digit_t *a, const digit_t *b, digit_t *c);
extern void fpsub747_asm(const digit_t *a, const digit_t *b, digit_t *c);

// Modular negation, a = -a mod p747
extern void fpneg747(digit_t *a);

// Modular division by two, c = a/2 mod p747.
void fpdiv2_747(const digit_t *a, digit_t *c);

// Modular correction to reduce field element a in [0, 2*p747-1] to [0, p747-1].
void fpcorrection747(digit_t *a);

// 747-bit Montgomery reduction, c = a mod p
void rdc_mont(const digit_t *a, digit_t *c);

// Field multiplication using Montgomery arithmetic, c = a*b*R^-1 mod p747, where R=2^768
void fpmul747_mont(const felm_t a, const felm_t b, felm_t c);
void mul747_asm(const felm_t a, const felm_t b, dfelm_t c);
void rdc747_asm(const dfelm_t ma, dfelm_t mc);

// Field squaring using Montgomery arithmetic, c = a*b*R^-1 mod p747, where R=2^768
void fpsqr747_mont(const felm_t ma, felm_t mc);

// Conversion to Montgomery representation
void to_mont(const felm_t a, felm_t mc);

// Conversion from Montgomery representation to standard representation
void from_mont(const felm_t ma, felm_t c);

// Field inversion, a = a^-1 in GF(p747)
void fpinv747_mont(felm_t a);

// Field inversion, a = a^-1 in GF(p747) using the binary GCD
void fpinv747_mont_bingcd(felm_t a);

// Chain to compute (p747-3)/4 using Montgomery arithmetic
void fpinv747_chain_mont(felm_t a);

/************ GF(p^2) arithmetic functions *************/

// Copy of a GF(p747^2) element, c = a
void fp2copy747(const f2elm_t a, f2elm_t c);

// Zeroing a GF(p747^2) element, a = 0
void fp2zero747(f2elm_t a);

// GF(p747^2) negation, a = -a in GF(p747^2)
void fp2neg747(f2elm_t a);

// GF(p747^2) addition, c = a+b in GF(p747^2)
extern void fp2add747(const f2elm_t a, const f2elm_t b, f2elm_t c);

// GF(p747^2) subtraction, c = a-b in GF(p747^2)
extern void fp2sub747(const f2elm_t a, const f2elm_t b, f2elm_t c);

// GF(p747^2) division by two, c = a/2  in GF(p747^2)
void fp2div2_747(const f2elm_t a, f2elm_t c);

// Modular correction, a = a in GF(p747^2)
void fp2correction747(f2elm_t a);

// GF(p747^2) squaring using Montgomery arithmetic, c = a^2 in GF(p747^2)
void fp2sqr747_mont(const f2elm_t a, f2elm_t c);

// GF(p747^2) multiplication using Montgomery arithmetic, c = a*b in GF(p747^2)
void fp2mul747_mont(const f2elm_t a, const f2elm_t b, f2elm_t c);

// Conversion of a GF(p747^2) element to Montgomery representation
void to_fp2mont(const f2elm_t a, f2elm_t mc);

// Conversion of a GF(p747^2) element from Montgomery representation to standard representation
void from_fp2mont(const f2elm_t ma, f2elm_t c);

// GF(p747^2) inversion using Montgomery arithmetic, a = (a0-i*a1)/(a0^2+a1^2)
void fp2inv747_mont(f2elm_t a);

// GF(p747^2) inversion, a = (a0-i*a1)/(a0^2+a1^2), GF(p747) inversion done using the binary GCD
void fp2inv747_mont_bingcd(f2elm_t a);

// n-way Montgomery inversion
void mont_n_way_inv(const f2elm_t *vec, const int n, f2elm_t *out);

/************ Elliptic curve and isogeny functions *************/

// Computes the j-invariant of a Montgomery curve with projective constant.
void j_inv(const f2elm_t A, const f2elm_t C, f2elm_t jinv);

// Simultaneous doubling and differential addition.
void xDBLADD(point_proj_t P, point_proj_t Q, const f2elm_t xPQ, const f2elm_t A24);

// Doubling of a Montgomery point in projective coordinates (X:Z).
void xDBL(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24);

// Computes [2^e](X:Z) on Montgomery curve with projective constant via e repeated doublings.
void xDBLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24, const int e);

// Differential addition.
void xADD(point_proj_t P, const point_proj_t Q, const f2elm_t xPQ);

// Computes the corresponding 4-isogeny of a projective Montgomery point (X4:Z4) of order 4.
void get_4_isog(const point_proj_t P, f2elm_t A24plus, f2elm_t C24, f2elm_t *coeff);

// Evaluates the isogeny at the point (X:Z) in the domain of the isogeny.
void eval_4_isog(point_proj_t P, f2elm_t *coeff);

// Tripling of a Montgomery point in projective coordinates (X:Z).
void xTPL(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24plus);

// Computes [3^e](X:Z) on Montgomery curve with projective constant via e repeated triplings.
void xTPLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24plus, const int e);

// Computes the corresponding 3-isogeny of a projective Montgomery point (X3:Z3) of order 3.
void get_3_isog(const point_proj_t P, f2elm_t A24minus, f2elm_t A24plus, f2elm_t *coeff);

// Computes the 3-isogeny R=phi(X:Z), given projective point (X3:Z3) of order 3 on a Montgomery curve and a point P with coefficients given in coeff.
void eval_3_isog(point_proj_t Q, const f2elm_t *coeff);

// 3-way simultaneous inversion
void inv_3_way(f2elm_t z1, f2elm_t z2, f2elm_t z3);

// Given the x-coordinates of P, Q, and R, returns the value A corresponding to the Montgomery curve E_A: y^2=x^3+A*x^2+x such that R=Q-P on E_A.
void get_A(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xR, f2elm_t A);

// 6-way simultaneous inversion
void inv_6_way(f2elm_t z1, f2elm_t z2, f2elm_t z3, f2elm_t z4, f2elm_t z5, f2elm_t z6);

// Criss cross operation for computing 5-isogenies
void criss_cross(f2elm_t alpha, f2elm_t beta, f2elm_t gamma, f2elm_t delta);

// Computes the 5-isogeny R=phi(X:Z), given two projective points P=(X3:Z3) and Pdbl=(X'3:Z'3) of order 3 on a Montgomery curve where Pdbl = [2]P
void eval_5_isog(const point_proj_t P, const point_proj_t Pdbl, point_proj_t R);

// Compute Montgomery curve projective coefficient from a projective point alpha of order 2 on the curve
void get_a_from_alpha(const point_proj_t alpha, f2elm_t A24plus, f2elm_t C24);

// Given the P, Q, and R, returns the value A24plus and C24 corresponding to the Montgomery curve E_A: Cy^2=Cx^3+A*x^2+Cx such that R=Q-P on E_A.
void get_A_projective(const point_proj_t P, const point_proj_t Q, const point_proj_t R, f2elm_t A24plus, f2elm_t C24);

// Simultaneous doubling and differential addition using projective curve coefficients.
void xDBLADD_AC24(point_proj_t P, point_proj_t Q, const point_proj_t xPQ, const f2elm_t A24plus, const f2elm_t C24);

// Quintupling of a Montgomery point in projective coordinates (X:Z).
void xQNTPL(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24);

// Computes [5^e](X:Z) on Montgomery curve with projective constant via e repeated quintuplings.
void xQNTPLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24, const int e);

void LADDER3PT(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xPQ, const digit_t *m, const unsigned int AliceOrBobOrEve, point_proj_t R, const f2elm_t A);

#endif
