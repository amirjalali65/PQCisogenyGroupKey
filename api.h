/********************************************************************************************
* Supersingular Isogeny Group Key Agreement Library
*
* Abstract: API header file for P747
* This file is the modified version of the P751 header file in SIKE implementation by Microsoft Research
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
// SECURITY NOTE: This library supports ephemeral group key exchange between three parties. It is NOT secure to use it with static keys.
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
// Output: the public key PublicKeyA consisting of 6 GF(p747^2) elements encoded in 1134 bytes.
int EphemeralKeyGeneration_A(const unsigned char* PrivateKeyA, unsigned char* PublicKeyA);

// Bob's ephemeral key-pair generation
// It produces a private key PrivateKeyB and computes the public key PublicKeyB.
// The private key is an integer in the range [0, 2^Floor(Log(2,3^153)) - 1], stored in 48 bytes.  
// The public key consists of 6 GF(p747^2) elements encoded in 1134 bytes.
int EphemeralKeyGeneration_B(const unsigned char* PrivateKeyB, unsigned char* PublicKeyB);

// Eve's ephemeral key-pair generation
// It produces a private key PrivateKeyC and computes the public key PublicKeyC.
// The private key is an integer in the range [0, 2^Floor(Log(2,5^105)) - 1], stored in 48 bytes.  
// The public key consists of 6 GF(p747^2) elements encoded in 1134 bytes.
int EphemeralKeyGeneration_C(const unsigned char* PrivateKeyC, unsigned char* PublicKeyC);

// Bob's shared public key generation 
// It produces a shared public key constructed between Alice and Bob using Bob's private key and Alice's Public key
// The private key is an integer in the range [0, 2^Floor(Log(2,3^153)) - 1], stored in 48 bytes.  
// The shared public key consists of 3 GF(p747^2) elements encoded in 567 bytes.
int BSharedPublicFromA(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, unsigned char* SharedPublicAB);

// Eve's shared secret key generation 
// It produces a shared public key constructed between Eve and Bob using Eve's private key and Bob's Public key
// It also generates the shared secret key from Eve's point of view using Alice and Bob shared public key
// The private key is an integer in the range [0, 2^Floor(Log(2,5^105)) - 1], stored in 48 bytes.  
// The shared public key consists of 3 GF(p747^2) elements encoded in 567 bytes.
// The shared secret key consists of one GF(p747^2) element encoded in 190 bytes.
int CSharedSecretFromB(const unsigned char* PrivateKeyC, const unsigned char* PublicKeyB, const unsigned char* SharedPublicAB, unsigned char* SharedPublicBC, unsigned char* SharedSecret);

// Alice's shared secret key generation 
// It produces a shared public key constructed between Alice and Eve using Alice's private key and Eve's Public key
// It also generates the shared secret key from Alice's point of view using Bob and Eve shared public key
// The private key is an integer in the range [0, 2^260 - 1], stored in 48 bytes.  
// The shared public key consists of 3 GF(p747^2) elements encoded in 567 bytes.
// The shared secret key consists of one GF(p747^2) element encoded in 190 bytes.
int ASharedSecretFromC(const unsigned char* PrivateKeyA, const unsigned char* PublicKeyC, const unsigned char* SharedPublicBC, unsigned char* SharedPublicAC, unsigned char* SharedSecret);

// Bob's shared secret key generation 
// It generates the shared secret key from Bob's point of view using Alice and Eve shared public key
// The private key is an integer in the range [0, 2^Floor(Log(2,3^153)) - 1], stored in 48 bytes.  
// The shared secret key consists of one GF(p747^2) element encoded in 190 bytes.
int BSharedSecretFromA(const unsigned char* PrivateKeyB, const unsigned char* SharedPublicAC, unsigned char* SharedSecret);

#endif
