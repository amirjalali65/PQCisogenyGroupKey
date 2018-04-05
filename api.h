/********************************************************************************************
* Supersingular Isogeny Group Key Agreement Library
*
* Abstract: API header file for P747
* This file is modified version of the P751 header file in SIKE implementation by Microsoft Research
* The modifications are based on the different prime and different funtions inside Group Key Library.
*
*   Modified by Amir Jalali             ajalali2016@fau.edu
*********************************************************************************************/  

#ifndef __P747_API_H__
#define __P747_API_H__

#include "config.h"
    

/*********************** Group Key Agreement API ***********************/

#define CRYPTO_SECRETKEYBYTES	         48    
#define CRYPTO_PUBLICKEYBYTES	         1134
#define CRYPTO_BYTES	                 190
#define CRYPTO_SHAREDPUBLICBYTES	     567    

// Algorithm name
#define CRYPTO_ALGNAME "SIGKp747"  

/*********************** Ephemeral group key exchange API ***********************/
/*
#define SIDH_SECRETKEYBYTES      48
#define SIDH_PUBLICKEYBYTES     564
#define SIDH_BYTES              188 
*/
// SECURITY NOTE: SIDH supports ephemeral Diffie-Hellman key exchange. It is NOT secure to use it with static keys.
// See "On the Security of Supersingular Isogeny Cryptosystems", S.D. Galbraith, C. Petit, B. Shani and Y.B. Ti, in ASIACRYPT 2016, 2016.
// Extended version available at: http://eprint.iacr.org/2016/859     

// Generation of Alice's secret key 
// Outputs random value in [0, 2^260 - 1] to be used as Alice's private key
void random_mod_order_A(unsigned char* random_digits);

// Generation of Bob's secret key 
// Outputs random value in [0, 2^Floor(Log(2,3^153)) - 1] to be used as Bob's private key
void random_mod_order_B(unsigned char* random_digits);

// Generation of Eve's secret key 
// Outputs random value in [0, 2^Floor(Log(2,5^105)) - 1] to be used as Eve's private key
void random_mod_order_C(unsigned char* random_digits);

// Alice's ephemeral public key generation
// Input:  a private key PrivateKeyA in the range [0, 2^260 - 1], stored in 47 bytes. 
// Output: the public key PublicKeyA consisting of 3 GF(p747^2) elements encoded in 564 bytes.
int EphemeralKeyGeneration_A(const unsigned char* PrivateKeyA, unsigned char* PublicKeyA);

// Bob's ephemeral key-pair generation
// It produces a private key PrivateKeyB and computes the public key PublicKeyB.
// The private key is an integer in the range [0, 2^Floor(Log(2,3^153)) - 1], stored in 48 bytes.  
// The public key consists of 3 GF(p747^2) elements encoded in 564 bytes.
int EphemeralKeyGeneration_B(const unsigned char* PrivateKeyB, unsigned char* PublicKeyB);

// Eve's ephemeral key-pair generation
// It produces a private key PrivateKeyC and computes the public key PublicKeyC.
// The private key is an integer in the range [0, 2^Floor(Log(2,5^105)) - 1], stored in 48 bytes.  
// The public key consists of 3 GF(p747^2) elements encoded in 564 bytes.
int EphemeralKeyGeneration_C(const unsigned char* PrivateKeyC, unsigned char* PublicKeyC);

int BSharedPublicFromA(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, unsigned char* SharedPublicAB);

int CSharedSecretFromB(const unsigned char* PrivateKeyC, const unsigned char* PublicKeyB, const unsigned char* SharedPublicAB, unsigned char* SharedPublicBC, unsigned char* SharedSecret);

int ASharedSecretFromC(const unsigned char* PrivateKeyA, const unsigned char* PublicKeyC, const unsigned char* SharedPublicBC, unsigned char* SharedPublicAC, unsigned char* SharedSecret);

int BSharedSecretFromA(const unsigned char* PrivateKeyB, const unsigned char* SharedPublicAC, unsigned char* SharedSecret);

#endif
