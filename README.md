# PQCisogenyGroupKey
The efficient implementation of supersingular isogeny 3-party group key agreement in C. 

This repository includes the ephemeral three-party key establishment from isogenies of supersingular elliptic curves which are assumed to be resistant against quantum attacks. 

This library is implemented based on the arithmetic implementation of SIDH library developed by [Microsoft Research](https://github.com/Microsoft/PQCrypto-SIDH). A set of new functionalities are implemented to compute and evaluate the isogenies of degree-5. Moreover, the protocol design and key exchange operations are completely different from 2-party key exchange.

## Finite Field
The proposed finite field is constructed over 747-bit prime providing 81-bit quantum security level. The proposed prime has a special form which makes it isogeny-friendly and Montgomery-friendly, resulting in efficient implementation of field and group operations. 

## API and Key Exchange Procedure

The proposed three-party key exchange is completed in after four passes. Suppose **A**, **B**, and **C** want to compute a shared group secret:

**A**, **B**, and **C** generate their secret keys using `random_mod_order_A()`, `random_mod_order_B()`, and `random_mod_order_C()`, respectively.

**A**, **B**, and **C** generate their public keys using `EphemeralKeyGeneration_A()`, `EphemeralKeyGeneration_B()`, and `EphemeralKeyGeneration_C()`, respectively.

1- **A** sends her public key to **B**. 

2- **B** computes a shared public key "Shared_AB" using his secret key and **A**'s public key by calling `BSharedPublicFromA()`. **B** sends his public key and the generated "Shared_AB" to **C**.

3- **C** computes a shared public key "Shared_BC" from  **B**'s public key along with the 3-party shared secret from "Shared_AB" using `CSharedSecretFromB()`. At this point, **C**'s computation is done and the shared secret is computed for this party. **C** send the generated "Shared_BC" and her public key to **A**.

4- **A** computes the 3-party shared secret using the "Shared_BC" received from **C**, and computes a shared public key "Shared_AC" using **C**'s public key using `ASharedSecretFromC()`. At this point, **A**'s computation is done and the shared secret is computed for this party. **A** sends the generated "Shared_AC" to **B** as the final pass. 
**B** computes the 3-party shared secret using "Shared_AC" by calling `BSharedSecretFromA()`. 

## Building Binary
This version of the library is implemented in C and supports different platforms. Simply use `make` in the terminal:
```sh
$ make 
```
## Contributors
Amir Jalali (ajalali[at]linkedin.com)
Reza Azarderakhsh (razarderakhsh@fau.edu)

