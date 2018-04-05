/********************************************************************************************
* Supersingular Isogeny Group Key Agreement Library
*
* Abstract: elliptic curve and isogeny functions
* This file is the modified version of the SIKE ec_isogeny file developed by Microsoft Research
* The added functions compute the isogeny of degree-5
*
*   Modified by Amir Jalali             ajalali2016@fau.edu
*********************************************************************************************/

#include "P747_internal.h"

void xDBL(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24)
{ // Doubling of a Montgomery point in projective coordinates (X:Z).
  // Input: projective Montgomery x-coordinates P = (X1:Z1), where x1=X1/Z1 and Montgomery curve constants A+2C and 4C.
  // Output: projective Montgomery x-coordinates Q = 2*P = (X2:Z2).
    f2elm_t t0, t1;
    
    fp2sub(P->X, P->Z, t0);                         // t0 = X1-Z1
    fp2add(P->X, P->Z, t1);                         // t1 = X1+Z1
    fp2sqr_mont(t0, t0);                            // t0 = (X1-Z1)^2 
    fp2sqr_mont(t1, t1);                            // t1 = (X1+Z1)^2 
    fp2mul_mont(C24, t0, Q->Z);                     // Z2 = C24*(X1-Z1)^2   
    fp2mul_mont(t1, Q->Z, Q->X);                    // X2 = C24*(X1-Z1)^2*(X1+Z1)^2
    fp2sub(t1, t0, t1);                             // t1 = (X1+Z1)^2-(X1-Z1)^2 
    fp2mul_mont(A24plus, t1, t0);                   // t0 = A24plus*[(X1+Z1)^2-(X1-Z1)^2]
    fp2add(Q->Z, t0, Q->Z);                         // Z2 = A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2
    fp2mul_mont(Q->Z, t1, Q->Z);                    // Z2 = [A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2]*[(X1+Z1)^2-(X1-Z1)^2]
}


void xDBLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24, const int e)
{ // Computes [2^e](X:Z) on Montgomery curve with projective constant via e repeated doublings.
  // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A+2C and 4C.
  // Output: projective Montgomery x-coordinates Q <- (2^e)*P.
    int i;
    
    copy_words((digit_t*)P, (digit_t*)Q, 2*2*NWORDS_FIELD);

    for (i = 0; i < e; i++) {
        xDBL(Q, Q, A24plus, C24);
    }
}


void get_4_isog(const point_proj_t P, f2elm_t A24plus, f2elm_t C24, f2elm_t* coeff)
{ // Computes the corresponding 4-isogeny of a projective Montgomery point (X4:Z4) of order 4.
  // Input:  projective point of order four P = (X4:Z4).
  // Output: the 4-isogenous Montgomery curve with projective coefficients A+2C/4C and the 3 coefficients 
  //         that are used to evaluate the isogeny at a point in eval_4_isog().
    
    fp2sub(P->X, P->Z, coeff[1]);                   // coeff[1] = X4-Z4
    fp2add(P->X, P->Z, coeff[2]);                   // coeff[2] = X4+Z4
    fp2sqr_mont(P->Z, coeff[0]);                    // coeff[0] = Z4^2
    fp2add(coeff[0], coeff[0], coeff[0]);           // coeff[0] = 2*Z4^2
    fp2sqr_mont(coeff[0], C24);                     // C24 = 4*Z4^4
    fp2add(coeff[0], coeff[0], coeff[0]);           // coeff[0] = 4*Z4^2
    fp2sqr_mont(P->X, A24plus);                     // A24plus = X4^2
    fp2add(A24plus, A24plus, A24plus);              // A24plus = 2*X4^2
    fp2sqr_mont(A24plus, A24plus);                  // A24plus = 4*X4^4
}


void eval_4_isog(point_proj_t P, f2elm_t* coeff)
{ // Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 4-isogeny phi defined 
  // by the 3 coefficients in coeff (computed in the function get_4_isog()).
  // Inputs: the coefficients defining the isogeny, and the projective point P = (X:Z).
  // Output: the projective point P = phi(P) = (X:Z) in the codomain. 
    f2elm_t t0, t1;
    
    fp2add(P->X, P->Z, t0);                         // t0 = X+Z
    fp2sub(P->X, P->Z, t1);                         // t1 = X-Z
    fp2mul_mont(t0, coeff[1], P->X);                // X = (X+Z)*coeff[1]
    fp2mul_mont(t1, coeff[2], P->Z);                // Z = (X-Z)*coeff[2]
    fp2mul_mont(t0, t1, t0);                        // t0 = (X+Z)*(X-Z)
    fp2mul_mont(t0, coeff[0], t0);                  // t0 = coeff[0]*(X+Z)*(X-Z)
    fp2add(P->X, P->Z, t1);                         // t1 = (X-Z)*coeff[2] + (X+Z)*coeff[1]
    fp2sub(P->X, P->Z, P->Z);                       // Z = (X-Z)*coeff[2] - (X+Z)*coeff[1]
    fp2sqr_mont(t1, t1);                            // t1 = [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
    fp2sqr_mont(P->Z, P->Z);                        // Z = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2
    fp2add(t1, t0, P->X);                           // X = coeff[0]*(X+Z)*(X-Z) + [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
    fp2sub(P->Z, t0, t0);                           // t0 = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2 - coeff[0]*(X+Z)*(X-Z)
    fp2mul_mont(P->X, t1, P->X);                    // Xfinal
    fp2mul_mont(P->Z, t0, P->Z);                    // Zfinal
}


void xTPL(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24plus)              
{ // Tripling of a Montgomery point in projective coordinates (X:Z).
  // Input: projective Montgomery x-coordinates P = (X:Z), where x=X/Z and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
  // Output: projective Montgomery x-coordinates Q = 3*P = (X3:Z3).
    f2elm_t t0, t1, t2, t3, t4, t5, t6;
                                    
    fp2sub(P->X, P->Z, t0);                         // t0 = X-Z 
    fp2sqr_mont(t0, t2);                            // t2 = (X-Z)^2           
    fp2add(P->X, P->Z, t1);                         // t1 = X+Z 
    fp2sqr_mont(t1, t3);                            // t3 = (X+Z)^2
    fp2add(t0, t1, t4);                             // t4 = 2*X
    fp2sub(t1, t0, t0);                             // t0 = 2*Z 
    fp2sqr_mont(t4, t1);                            // t1 = 4*X^2
    fp2sub(t1, t3, t1);                             // t1 = 4*X^2 - (X+Z)^2 
    fp2sub(t1, t2, t1);                             // t1 = 4*X^2 - (X+Z)^2 - (X-Z)^2
    fp2mul_mont(t3, A24plus, t5);                   // t5 = A24plus*(X+Z)^2 
    fp2mul_mont(t3, t5, t3);                        // t3 = A24plus*(X+Z)^3
    fp2mul_mont(A24minus, t2, t6);                  // t6 = A24minus*(X-Z)^2
    fp2mul_mont(t2, t6, t2);                        // t2 = A24minus*(X-Z)^3
    fp2sub(t2, t3, t3);                             // t3 = A24minus*(X-Z)^3 - coeff*(X+Z)^3
    fp2sub(t5, t6, t2);                             // t2 = A24plus*(X+Z)^2 - A24minus*(X-Z)^2
    fp2mul_mont(t1, t2, t1);                        // t1 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
    fp2add(t3, t1, t2);                             // t2 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2] + A24minus*(X-Z)^3 - coeff*(X+Z)^3
    fp2sqr_mont(t2, t2);                            // t2 = t2^2
    fp2mul_mont(t4, t2, Q->X);                      // X3 = 2*X*t2
    fp2sub(t3, t1, t1);                             // t1 = A24minus*(X-Z)^3 - A24plus*(X+Z)^3 - [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
    fp2sqr_mont(t1, t1);                            // t1 = t1^2
    fp2mul_mont(t0, t1, Q->Z);                      // Z3 = 2*Z*t1
}


void xTPLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24plus, const int e)
{ // Computes [3^e](X:Z) on Montgomery curve with projective constant via e repeated triplings.
  // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
  // Output: projective Montgomery x-coordinates Q <- (3^e)*P.
    int i;
        
    copy_words((digit_t*)P, (digit_t*)Q, 2*2*NWORDS_FIELD);

    for (i = 0; i < e; i++) {
        xTPL(Q, Q, A24minus, A24plus);
    }
}


void get_3_isog(const point_proj_t P, f2elm_t A24minus, f2elm_t A24plus, f2elm_t* coeff)
{ // Computes the corresponding 3-isogeny of a projective Montgomery point (X3:Z3) of order 3.
  // Input:  projective point of order three P = (X3:Z3).
  // Output: the 3-isogenous Montgomery curve with projective coefficient A/C. 
    f2elm_t t0, t1, t2, t3, t4;
    
    fp2sub(P->X, P->Z, coeff[0]);                   // coeff0 = X-Z
    fp2sqr_mont(coeff[0], t0);                      // t0 = (X-Z)^2
    fp2add(P->X, P->Z, coeff[1]);                   // coeff1 = X+Z
    fp2sqr_mont(coeff[1], t1);                      // t1 = (X+Z)^2
    fp2add(t0, t1, t2);                             // t2 = (X+Z)^2 + (X-Z)^2
    fp2add(coeff[0], coeff[1], t3);                 // t3 = 2*X
    fp2sqr_mont(t3, t3);                            // t3 = 4*X^2
    fp2sub(t3, t2, t3);                             // t3 = 4*X^2 - (X+Z)^2 - (X-Z)^2 
    fp2add(t1, t3, t2);                             // t2 = 4*X^2 - (X-Z)^2 
    fp2add(t3, t0, t3);                             // t3 = 4*X^2 - (X+Z)^2
    fp2add(t0, t3, t4);                             // t4 = 4*X^2 - (X+Z)^2 + (X-Z)^2 
    fp2add(t4, t4, t4);                             // t4 = 2(4*X^2 - (X+Z)^2 + (X-Z)^2) 
    fp2add(t1, t4, t4);                             // t4 = 8*X^2 - (X+Z)^2 + 2*(X-Z)^2
    fp2mul_mont(t2, t4, A24minus);                  // A24minus = [4*X^2 - (X-Z)^2]*[8*X^2 - (X+Z)^2 + 2*(X-Z)^2]
    fp2add(t1, t2, t4);                             // t4 = 4*X^2 + (X+Z)^2 - (X-Z)^2
    fp2add(t4, t4, t4);                             // t4 = 2(4*X^2 + (X+Z)^2 - (X-Z)^2) 
    fp2add(t0, t4, t4);                             // t4 = 8*X^2 + 2*(X+Z)^2 - (X-Z)^2
    fp2mul_mont(t3, t4, t4);                        // t4 = [4*X^2 - (X+Z)^2]*[8*X^2 + 2*(X+Z)^2 - (X-Z)^2]
    fp2sub(t4, A24minus, t0);                       // t0 = [4*X^2 - (X+Z)^2]*[8*X^2 + 2*(X+Z)^2 - (X-Z)^2] - [4*X^2 - (X-Z)^2]*[8*X^2 - (X+Z)^2 + 2*(X-Z)^2] 
    fp2add(A24minus, t0, A24plus);                  // A24plus = 8*X^2 - (X+Z)^2 + 2*(X-Z)^2
}


void eval_3_isog(point_proj_t Q, const f2elm_t* coeff)
{ // Computes the 3-isogeny R=phi(X:Z), given projective point (X3:Z3) of order 3 on a Montgomery curve and 
  // a point P with 2 coefficients in coeff (computed in the function get_3_isog()).
  // Inputs: projective points P = (X3:Z3) and Q = (X:Z).
  // Output: the projective point Q <- phi(Q) = (X3:Z3). 
    f2elm_t t0, t1, t2;

    fp2add(Q->X, Q->Z, t0);                       // t0 = X+Z
    fp2sub(Q->X, Q->Z, t1);                       // t1 = X-Z
    fp2mul_mont(t0, coeff[0], t0);                // t0 = coeff0*(X+Z)
    fp2mul_mont(t1, coeff[1], t1);                // t1 = coeff1*(X-Z)
    fp2add(t0, t1, t2);                           // t2 = coeff0*(X-Z) + coeff1*(X+Z)
    fp2sub(t1, t0, t0);                           // t0 = coeff0*(X-Z) - coeff1*(X+Z)
    fp2sqr_mont(t2, t2);                          // t2 = [coeff0*(X-Z) + coeff1*(X+Z)]^2
    fp2sqr_mont(t0, t0);                          // t1 = [coeff0*(X-Z) - coeff1*(X+Z)]^2
    fp2mul_mont(Q->X, t2, Q->X);                  // X3final = X*[coeff0*(X-Z) + coeff1*(X+Z)]^2        
    fp2mul_mont(Q->Z, t0, Q->Z);                  // Z3final = Z*[coeff0*(X-Z) - coeff1*(X+Z)]^2
}


void inv_3_way(f2elm_t z1, f2elm_t z2, f2elm_t z3)
{ // 3-way simultaneous inversion
  // Input:  z1,z2,z3
  // Output: 1/z1,1/z2,1/z3 (override inputs).
    f2elm_t t0, t1, t2, t3;

    fp2mul_mont(z1, z2, t0);                      // t0 = z1*z2
    fp2mul_mont(z3, t0, t1);                      // t1 = z1*z2*z3
    fp2inv_mont(t1);                              // t1 = 1/(z1*z2*z3)
    fp2mul_mont(z3, t1, t2);                      // t2 = 1/(z1*z2) 
    fp2mul_mont(t2, z2, t3);                      // t3 = 1/z1
    fp2mul_mont(t2, z1, z2);                      // z2 = 1/z2
    fp2mul_mont(t0, t1, z3);                      // z3 = 1/z3
    fp2copy(t3, z1);                              // z1 = 1/z1
}

void inv_6_way(f2elm_t z1, f2elm_t z2, f2elm_t z3, f2elm_t z4, f2elm_t z5, f2elm_t z6)
{ // 3-way simultaneous inversion
  // Input:  z1,z2,z3, z4, z5, z6
  // Output: 1/z1,1/z2,1/z3,1/z4,1/z5,1/z6 (override inputs).
    f2elm_t t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    fp2mul_mont(z1, z2, t0);                      // t0 = z1*z2
    fp2mul_mont(z3, t0, t1);                      // t1 = z1*z2*z3
    fp2mul_mont(z4, t1, t2);                      // t2 = z1*z2*z3*z4
    fp2mul_mont(z5, t2, t3);                      // t3 = z1*z2*z3*z4*z5
    fp2mul_mont(z6, t3, t4);                      // t4 = z1*z2*z3*z4*z5*z6
    fp2inv_mont(t4);                           // t4 = 1/(z1*z2*z3*z4*z5*z6)
    fp2mul_mont(t4, t3, t5);                      // t5 = 1/z6
    fp2mul_mont(t4, z6, t4);                      // t4 = 1/(z1*z2*z3*z4*z5)
    fp2mul_mont(t4, t2, t6);                      // t6 = 1/z5
    fp2mul_mont(t4, z5, t4);                      // t4 = 1/(z1*z2*z3*z4)
    fp2mul_mont(t4, t1, t8);                      // t8 = 1/z4
    fp2mul_mont(t4, z4, t4);                      // t4 = 1/(z1*z2*z3)
    fp2mul_mont(t4, t0, t9);                      // t9 = 1/z3
    fp2mul_mont(t4, z3, t4);                      // t4 = 1/(z1*z2)
    fp2mul_mont(t4, z2, t10);                     // t10 = 1/z1
    fp2mul_mont(t4, z1, z2);                      // z2 = 1/z2
    fp2copy(t10, z1);                             // z1 = 1/z1  
    fp2copy(t9, z3);                              // z3 = 1/z3
    fp2copy(t8, z4);                              // z4 = 1/z4
    fp2copy(t6, z5);                              // z5 = 1/z5
    fp2copy(t5, z6);                              // z6 = 1/z6

}

void get_A(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xR, f2elm_t A)
{ // Given the x-coordinates of P, Q, and R, returns the value A corresponding to the Montgomery curve E_A: y^2=x^3+A*x^2+x such that R=Q-P on E_A.
  // Input:  the x-coordinates xP, xQ, and xR of the points P, Q and R.
  // Output: the coefficient A corresponding to the curve E_A: y^2=x^3+A*x^2+x.
    f2elm_t t0, t1, one = {0};
    
    fpcopy((digit_t*)&Montgomery_one, one[0]);
    fp2add(xP, xQ, t1);                           // t1 = xP+xQ
    fp2mul_mont(xP, xQ, t0);                      // t0 = xP*xQ
    fp2mul_mont(xR, t1, A);                       // A = xR*t1
    fp2add(t0, A, A);                             // A = A+t0
    fp2mul_mont(t0, xR, t0);                      // t0 = t0*xR
    fp2sub(A, one, A);                            // A = A-1
    fp2add(t0, t0, t0);                           // t0 = t0+t0
    fp2add(t1, xR, t1);                           // t1 = t1+xR
    fp2add(t0, t0, t0);                           // t0 = t0+t0
    fp2sqr_mont(A, A);                            // A = A^2
    fp2inv_mont(t0);                              // t0 = 1/t0
    fp2mul_mont(A, t0, A);                        // A = A*t0
    fp2sub(A, t1, A);                             // Afinal = A-t1
}


void j_inv(const f2elm_t A, const f2elm_t C, f2elm_t jinv)
{ // Computes the j-invariant of a Montgomery curve with projective constant.
  // Input: A,C in GF(p^2).
  // Output: j=256*(A^2-3*C^2)^3/(C^4*(A^2-4*C^2)), which is the j-invariant of the Montgomery curve B*y^2=x^3+(A/C)*x^2+x or (equivalently) j-invariant of B'*y^2=C*x^3+A*x^2+C*x.
    f2elm_t t0, t1;
    
    fp2sqr_mont(A, jinv);                           // jinv = A^2        
    fp2sqr_mont(C, t1);                             // t1 = C^2
    fp2add(t1, t1, t0);                             // t0 = t1+t1
    fp2sub(jinv, t0, t0);                           // t0 = jinv-t0
    fp2sub(t0, t1, t0);                             // t0 = t0-t1
    fp2sub(t0, t1, jinv);                           // jinv = t0-t1
    fp2sqr_mont(t1, t1);                            // t1 = t1^2
    fp2mul_mont(jinv, t1, jinv);                    // jinv = jinv*t1
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2sqr_mont(t0, t1);                            // t1 = t0^2
    fp2mul_mont(t0, t1, t0);                        // t0 = t0*t1
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2inv_mont(jinv);                              // jinv = 1/jinv 
    fp2mul_mont(jinv, t0, jinv);                    // jinv = t0*jinv
}


void xDBLADD(point_proj_t P, point_proj_t Q, const f2elm_t xPQ, const f2elm_t A24)
{ // Simultaneous doubling and differential addition.
  // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.
  // Output: projective Montgomery points P <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P, and Q <- P+Q = (XQP:ZQP) such that = x(Q+P)=XQP/ZQP. 
    f2elm_t t0, t1, t2;

    fp2add(P->X, P->Z, t0);                         // t0 = XP+ZP
    fp2sub(P->X, P->Z, t1);                         // t1 = XP-ZP
    fp2sqr_mont(t0, P->X);                          // XP = (XP+ZP)^2
    fp2sub(Q->X, Q->Z, t2);                         // t2 = XQ-ZQ
    fp2correction(t2);
    fp2add(Q->X, Q->Z, Q->X);                       // XQ = XQ+ZQ
    fp2mul_mont(t0, t2, t0);                        // t0 = (XP+ZP)*(XQ-ZQ)
    fp2sqr_mont(t1, P->Z);                          // ZP = (XP-ZP)^2
    fp2mul_mont(t1, Q->X, t1);                      // t1 = (XP-ZP)*(XQ+ZQ)
    fp2sub(P->X, P->Z, t2);                         // t2 = (XP+ZP)^2-(XP-ZP)^2
    fp2mul_mont(P->X, P->Z, P->X);                  // XP = (XP+ZP)^2*(XP-ZP)^2
    fp2mul_mont(t2, A24, Q->X);                     // XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sub(t0, t1, Q->Z);                           // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
    fp2add(Q->X, P->Z, P->Z);                       // ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
    fp2add(t0, t1, Q->X);                           // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
    fp2mul_mont(P->Z, t2, P->Z);                    // ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sqr_mont(Q->Z, Q->Z);                        // ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
    fp2sqr_mont(Q->X, Q->X);                        // XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
    fp2mul_mont(Q->Z, xPQ, Q->Z);                   // ZQ = xPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
}

void xDBLADD_AC24(point_proj_t P, point_proj_t Q, const point_proj_t PQ, const f2elm_t A24plus, const f2elm_t C24)
{   // Compute Double and Add using Montgomery formulas considering curve projective coefficient
    // A24plus = (A + 2C) and C24 = 4C P <- 2P and Q <- P + Q
	f2elm_t t0, t1, t2, pz, px;

	fp2copy(PQ->X, px);
	fp2copy(PQ->Z, pz);
	fp2add(P->X, P->Z, t0);                         // t0 = XP+ZP
	fp2sub(P->X, P->Z, t1);                         // t1 = XP-ZP
	fp2sqr_mont(t0, P->X);                          // XP = (XP+ZP)^2
	fp2sub(Q->X, Q->Z, t2);                         // t2 = XQ-ZQ
	fp2correction(t2);
	fp2add(Q->X, Q->Z, Q->X);                       // XQ = XQ+ZQ
	fp2mul_mont(t0, t2, t0);                        // t0 = (XP+ZP)*(XQ-ZQ)
	fp2sqr_mont(t1, P->Z);                          // ZP = (XP-ZP)^2
	fp2mul_mont(t1, Q->X, t1);                      // t1 = (XP-ZP)*(XQ+ZQ)
	fp2sub(P->X, P->Z, t2);                         // t2 = (XP+ZP)^2-(XP-ZP)^2
	fp2mul_mont(P->X, P->Z, P->X);                  // XP = (XP+ZP)^2*(XP-ZP)^2
	fp2mul_mont(P->X, C24, P->X);					// XP = C24*(XP+ZP)^2*(XP-ZP)^2
	fp2mul_mont(t2, A24plus, Q->X);                 // XQ = A24plus*[(XP+ZP)^2-(XP-ZP)^2]
	fp2sub(t0, t1, Q->Z);                           // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
	fp2mul_mont(P->Z, C24, P->Z);					// ZP = C24*(XP-ZP)^2
	fp2add(Q->X, P->Z, P->Z);                       // ZP = A24plus*[(XP+ZP)^2-(XP-ZP)^2]+C24*(XP-ZP)^2
	fp2add(t0, t1, Q->X);                           // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
	fp2mul_mont(P->Z, t2, P->Z);                    // ZP = [A24plus*[(XP+ZP)^2-(XP-ZP)^2]+C24*(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]
	fp2sqr_mont(Q->Z, Q->Z);                        // ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
	fp2sqr_mont(Q->X, Q->X);                        // XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
	fp2mul_mont(Q->X, pz, Q->X);					// XQ = ZPQ*[(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
	fp2mul_mont(Q->Z, px, Q->Z);	                // ZQ = XPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
}
static void swap_points(point_proj_t P, point_proj_t Q, const digit_t option)
{ // Swap points.
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    digit_t temp;
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++) {
        temp = option & (P->X[0][i] ^ Q->X[0][i]);
        P->X[0][i] = temp ^ P->X[0][i]; 
        Q->X[0][i] = temp ^ Q->X[0][i]; 
        temp = option & (P->Z[0][i] ^ Q->Z[0][i]);
        P->Z[0][i] = temp ^ P->Z[0][i]; 
        Q->Z[0][i] = temp ^ Q->Z[0][i]; 
        temp = option & (P->X[1][i] ^ Q->X[1][i]);
        P->X[1][i] = temp ^ P->X[1][i]; 
        Q->X[1][i] = temp ^ Q->X[1][i]; 
        temp = option & (P->Z[1][i] ^ Q->Z[1][i]);
        P->Z[1][i] = temp ^ P->Z[1][i]; 
        Q->Z[1][i] = temp ^ Q->Z[1][i]; 
    }
}

void LADDER3PT(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xPQ, const digit_t* m, const unsigned int AliceOrBobOrEve, point_proj_t R, const f2elm_t A)
{
    point_proj_t R0 = {0}, R2 = {0};
    f2elm_t A24 = {0};
    digit_t mask;
    int i, nbits, bit, swap, prevbit = 0;

    if (AliceOrBobOrEve == ALICE) {
        nbits = OALICE_BITS;
    } else if (AliceOrBobOrEve == BOB) {
        nbits = OBOB_BITS;
    } else {
        nbits = OEVE_BITS;
    }

    // Initializing constant
    fpcopy((digit_t*)&Montgomery_one, A24[0]);
    fp2add(A24, A24, A24);
    fp2add(A, A24, A24);
    fp2div2(A24, A24);  
    fp2div2(A24, A24); // A24 = (A+2)/4

    // Initializing points
    fp2copy(xQ, R0->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R0->Z);
    fp2copy(xPQ, R2->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R2->Z);
    fp2copy(xP, R->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R->Z);
    fpzero((digit_t*)(R->Z)[1]);

    // Main loop
    for (i = 0; i < nbits; i++) {
        bit = (m[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(R, R2, mask);
        xDBLADD(R0, R2, R->X, A24);
        fp2mul_mont(R2->X, R->Z, R2->X);
    }
}

void xQNTPL_ladder(const point_proj_t P, point_proj_t R, const f2elm_t A24plus, const f2elm_t C24)
{
	point_proj_t Q, R0;
	
	fp2copy(P->X, R0->X);
    fp2copy(P->Z, R0->Z);
	xDBL(R0, Q, A24plus, C24);                       // Q = 2*P
	xDBLADD_AC24(R0, Q, P, A24plus, C24);		 // R0 = 2*P, Q = 3*P
	xDBLADD_AC24(R0, Q, P, A24plus, C24);		 // R0 = 4*P, Q = 5*P
    fp2copy(Q->X, R->X);
    fp2copy(Q->Z, R->Z);
}

void criss_cross(f2elm_t alpha, f2elm_t beta, f2elm_t gamma, f2elm_t delta)
{   // CrissCross operation corresponds to Costello and Hisil method
    // Input: alpha, beta, gamma, and delta in K
    // Output: alpha = (alpha*delta + beta*gamma)
    //         beta = (alpha*delta - beta*gamma)

    f2elm_t t0, t1;
    fp2mul_mont(alpha, delta, t0);
    fp2mul_mont(beta, gamma, t1);
    fp2add(t0, t1, alpha);
    fp2sub(t0, t1, beta);
}

void eval_5_isog(const point_proj_t P, const point_proj_t Pdbl, point_proj_t Q)
{   // Compute the image of a point Q over a 5-isogeny 
    // Inputs are two points P  coordinates, and a point Q
    // Output is a point Q <- phi(Q)
    f2elm_t X_hat, Z_hat, t0, t1, t2, t3, t4, t5;

	fp2add(P->X, P->Z, t2);
	fp2sub(P->X, P->Z, t3);
	fp2add(Pdbl->X, Pdbl->Z, t4);
	fp2sub(Pdbl->X, Pdbl->Z, t5);
	
	fp2add(Q->X, Q->Z, X_hat);
    fp2sub(Q->X, Q->Z, Z_hat);
    criss_cross(t2, t3, X_hat, Z_hat);
    criss_cross(t4, t5, X_hat, Z_hat);    
    fp2mul_mont(t4, t2, t0);
    fp2mul_mont(t5, t3, t1);
    fp2sqr_mont(t0, t0); 
    fp2sqr_mont(t1, t1);
    fp2mul_mont(t0, Q->X, Q->X);
    fp2mul_mont(t1, Q->Z, Q->Z);
}

void get_a_from_alpha(const point_proj_t alpha, f2elm_t A24plus, f2elm_t C24)
{   // Compute the curve coefficient A24plus = (A + 2C) and C24 = 4C using
    // 2-torsion point alpha on the curve E. This function is used to retrieve 
    // curve coefficients from the image of 2-torsion points on the new curve
    // Input: the projective coordinates of alpha
    // Output: Curve coefficients (A + 2C), 4C
	fp2sub(alpha->X, alpha->Z, A24plus);		// A24plus = (Xalpha - Zalpha)
	fp2sqr_mont(A24plus, A24plus);				// A24plus = (Xalpha - Zalpha)^2
	fp2correction(A24plus);
	fp2add(alpha->X, alpha->Z, C24);			// C24 = Xalpha + Zalpha
	fp2sqr_mont(C24, C24);						// C24 = (Xalpha + Zalpha)^2
	fp2sub(A24plus, C24, C24);					// C24 = (Xalpha - Zalpha)^2 - (Xalpha + Zalpha)^2
	fp2correction(C24);
}

void get_A_projective(const point_proj_t P, const point_proj_t Q, const point_proj_t R, f2elm_t A24plus, f2elm_t C24)
{   // Compute the curve coefficient A24plus = (A + 2C) and C24 = 4C using 
    // three projective points P, Q, and R where R = Q - P 
    // This is the projective version of get_A function 
	f2elm_t t0, t1, t2, t3, t4, t5, t6, t7, t8;
	fp2mul_mont(P->X, Q->X, t0);					// t0 = XP*XQ
	fp2mul_mont(t0, R->X, t0);						// t0 = XP*XQ*XR
	fp2mul_mont(P->Z, Q->Z, t1);					// t1 = ZP*ZQ
	fp2mul_mont(t1, R->Z, t1);						// t1 = ZP*ZQ*ZR
	fp2add(t0, t0, t0);								// t0 = 2*(XP*XQ*XR)
	fp2add(t1, t1, t8);								// t8 = 2*(ZP*ZQ*ZR)
	fp2mul_mont(t8, t0, C24);						// C24 = 4*(XP*XQ*XR)*(ZP*ZQ*ZR)
	fp2mul_mont(P->X, R->Z, t2);					// t2 = XP*ZR
	fp2mul_mont(Q->X, P->Z, t3);					// t3 = XQ*ZP
	fp2mul_mont(R->X, Q->Z, t4);					// t4 = XR*ZQ
	fp2mul_mont(t2, Q->Z, t5);						// t5 = XP*ZR*ZQ
	fp2mul_mont(t3, R->Z, t6);						// t6 = XQ*ZP*ZR
	fp2mul_mont(t4, P->Z, t7);						// t7 = XR*ZQ*ZP
	fp2add(t5, t6, t6);
	fp2add(t6, t7, t6);								// t6 = XP*ZR*ZQ + XQ*ZP*ZR + XR*ZQ*ZP
	fp2add(t0, t0, t0);								// t0 = 4*(XP*XQ*XR)
	fp2mul_mont(t0, t6, t6);						// t6 = 4*(XP*XQ*XR)*(XP*ZR*ZQ + XQ*ZP*ZR + XR*ZQ*ZP)
	fp2neg(t6);										// t6 = -4*(XP*XQ*XR)*(XP*ZR*ZQ + XQ*ZP*ZR + XR*ZQ*ZP)
	fp2mul_mont(t2, Q->X, t2);						// t2 = XP*XQ*ZR
	fp2mul_mont(t3, R->X, t3);						// t3 = XQ*XR*ZP
	fp2mul_mont(t4, P->X, t4);						// t4 = XR*XP*ZQ
	fp2add(t2, t3, t2);								// t2 = XP*XQ*ZR + XQ*XR*ZP
	fp2add(t2, t4, t2);								// t2 = XP*XQ*ZR + XQ*XR*ZP + XR*XP*ZQ
	fp2sub(t2, t1, t2);								// t2 = XP*XQ*ZR + XQ*XR*ZP + XR*XP*ZQ - ZP*ZQ*ZR
	fp2sqr_mont(t2, t2);							// t2 = (XP*XQ*ZR + XQ*XR*ZP + XR*XP*ZQ - ZP*ZQ*ZR)^2
	fp2add(t2, t6, t2);								// t2 = (XP*XQ*ZR + XQ*XR*ZP + XR*XP*ZQ - ZP*ZQ*ZR)^2 -4*(XP*XQ*XR)*(XP*ZR*ZQ + XQ*ZP*ZR + XR*ZQ*ZP)
	fp2add(C24, C24, t0);							// t0 = 8*(XP*XQ*XR)*(ZP*ZQ*ZR)
	fp2add(t2, t0, A24plus);						// A24plus = (XP*XQ*ZR + XQ*XR*ZP + XR*XP*ZQ - ZP*ZQ*ZR)^2 -4*(XP*XQ*XR)*(XP*ZR*ZQ + XQ*ZP*ZR + XR*ZQ*ZP) + 8*(XP*XQ*XR)*(ZP*ZQ*ZR)
	fp2add(t0, t0, C24);							// C24 = 16*(XP*XQ*XR)*(ZP*ZQ*ZR)
	fp2correction(A24plus);
	fp2correction(C24);
}

void xQNTPLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24, const int e)
{ // Computes [5^e](X:Z) on Montgomery curve with projective constant via e repeated triplings.
  // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
  // Output: projective Montgomery x-coordinates Q <- (5^e)*P.
    int i;
        
    copy_words((digit_t*)P, (digit_t*)Q, 2*2*NWORDS_FIELD);

    for (i = 0; i < e; i++) {
        xQNTPL_ladder(Q, Q, A24plus, C24);
    }
}
