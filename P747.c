/********************************************************************************************
* Supersingular Isogeny Group Key Agreement Library
*
* Abstract: supersingular isogeny parameters and generation of functions for P747
* The structure of this file is based on P751 file inside the SIKE library developed by Microsoft Research
*
* Modified by Amir Jalali           ajalali2016@fau.edu
*********************************************************************************************/  

#include "P747_internal.h"

// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 747-bit field element is represented with Ceil(747 / 64) = 12 64-bit digits or Ceil(747 / 32) = 24 32-bit digits.

// Curve isogeny system "SIGKp747". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p747^2), where A=0, B=1, C=1 and   
//
         
const uint64_t p747[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xA047C05513A78CEF, 0xB26D38D6F3B27CC3,
                                                     0x0B40709FDFCF993C, 0x7B53A41A27D10162, 0xA4DC87C4B86348CC, 0x5F33FCB0E1016AA2, 0xD42A27A9491431BC, 0x000004EE30756B8B };
const uint64_t p747p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xA047C05513A78CF0, 0xB26D38D6F3B27CC3,
                                                     0x0B40709FDFCF993C, 0x7B53A41A27D10162, 0xA4DC87C4B86348CC, 0x5F33FCB0E1016AA2, 0xD42A27A9491431BC, 0x000004EE30756B8B };
const uint64_t p747x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x408F80AA274F19DF, 0x64DA71ADE764F987, 
                                                     0x1680E13FBF9F3279, 0xF6A748344FA202C4, 0x49B90F8970C69198, 0xBE67F961C202D545, 0xA8544F5292286378, 0x000009DC60EAD717 }; 
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000010 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0x0000000000000000, 0x6A51808385F9D8A3, 0xD30435F904DF3586, 0x1A029D363E05965E, 0x0005A76A2991D4FE };
// Order of Eve's subgroup
const uint64_t Eve_order[NWORDS64_ORDER]         = { 0x0000000000000000, 0xFB20C48AD00585E5, 0x3ED67BA06FA1853C, 0x3E12F2967B66737E, 0x000DF3D5E9BC0F65 };
// Alice's generator values {XPA0 + XPA1*i, XQA0, XRA0 + XRA1*i} in GF(p747^2), expressed in normal representation
const uint64_t A_gen[5 * NWORDS64_FIELD]         = { 0x146A64BF56F93A7C, 0xD2834AEB7FAFAD64, 0xA813E25F64724ECA, 0x263CAEFDCFBC9279, 0x94D8C091FBE820C2, 0xF3FD5F9EB76FD467,
													 0x53FAD378BD2824EA, 0xCA4BF0D29F09B061, 0x3A3B1CC4F0B926F7, 0x768CC2152752FA5E, 0xED1D40B964662E78, 0x00000120A5B313BA,   // XPA0
													 0x1E6A90AEC79F4435, 0x636DCCE289A19199, 0x25A5C1A36709082C, 0xA1F0B1F01A226759, 0x810D8C4C978BD734, 0x175A804F0A2D4C37,
													 0x05956FCE365275A1, 0x4C0DEA39E9FA3121, 0xC09528C4A8DF299D, 0x8DC034AA3577B198, 0x60D67E17D7F8C860, 0x000000B9D6998639,   // XPA1
													 0xF328FA10F91C45F0, 0xE5A055346EA60C70, 0xDFDA473DEB9931C3, 0x4633D775F2407AC6, 0x3E21A2C1599493C4, 0xB24A13A85E621EE0,
													 0xCDEA5A68DCD0B2F2, 0xA6D518EDB17B32A4, 0xC7D196FA85A9E39D, 0x1331646D73439934, 0x310117A81F0143FA, 0x0000021D6762FF18,   // XQA0
													 0x124024C5480C2696, 0x7290343A5864802D, 0x3B7A746AE11871BA, 0xA3969F3C2099AA85, 0x5674927D92F1BCC1, 0x9FB3BCA6B6AC1ECC,
													 0xB11FDAF64CEF67EE, 0x64E250AC0B9FA8F2, 0x6CCDDDD25F56A1E6, 0xC2F7EFE77827FB7D, 0x9578C5F557EB62D9, 0x000004723AC260D5,   // XRA0
													 0x74DC8E0FD9052C39, 0x78A4DED7648B4B52, 0x19BD6A179F43E717, 0x821C4EAC5AFC0DAD, 0xF896042098451E78, 0xD3553C0D99F4933B,
													 0xA3BCC31111792301, 0x4F1AB67D511326EE, 0x54452EAD8482B25F, 0x1B99283D8D928DF4, 0x9003A7877DAE4AF9, 0x00000270E6E06619 }; // XRA1
// Bob's generator values {XPB0 + XPB1*i, XQB0, XRB0 + XRB1*i} in GF(p747^2), expressed in normal representation
const uint64_t B_gen[5 * NWORDS64_FIELD]         = { 0x9EE4AC530EA02812, 0x92C080440723255B, 0x662C55DBA078BBE3, 0x48B22316211DBAD4, 0xDE356317C914373B, 0xF78ED441F1DF05D0,
													 0x3111DFCCECCBD48C, 0x6720B43876BD4C8C, 0x99EE79475E08834F, 0x11DBD2F070A76299, 0x2F589404C5A6A8B2, 0x000004C377C95424,   // XPB0
													 0x1C7D4234E5FDCC74, 0x4DDAC3F7ADC53F78, 0xA84B1D9E5F46AB8E, 0xFC50A0657655B9C2, 0xF888E86F40EABDC1, 0xA496C18DA958AB38,
													 0x433E22772CD614FC, 0x4C2B0917B6D87723, 0xFB5E98C36C86388F, 0x18170BDA0CD711F6, 0x65A1BFA3BA76ADF9, 0x0000009A1D4C464D,   // XPB1
													 0x6D1BAEDCF00F6471, 0x448D26F2BD69042E, 0x35CE3DF10EF1B224, 0xD6CECADBA5451FEC, 0x268DBAFCFFB3499E, 0xBBEB5F0C9DEF37CA,
													 0x5B9F9109AA203E96, 0x65807C9E65B64504, 0xF302FE3DCF71BE79, 0x18073BC4322D75EB, 0xF606FD0C2F8FC5B1, 0x00000160C487D33C,   // XQB0
													 0xB22675E3A91F0902, 0xCDA1170DDD175E4F, 0x4DC79EFD82ECC131, 0x527554433D0294F5, 0x3EA091E8E417E852, 0xFFA76D7A98CDC144,
													 0x333A0B67E8B38716, 0x4DA35A16E089A0E6, 0xEAB4838DAD241FC4, 0x2BB1E64C0B454D30, 0xC3B2FB82628FA06E, 0x0000031431B95584,   // XRB0
													 0x364F7B32FAE86420, 0x4263E9F2477348EE, 0x2B81A33361D8687A, 0x64911A7CD8084228, 0x66AFB18A486140E8, 0xF2184390441F7512,
													 0xB5DE065CCD4F116E, 0xA43BDE0F0B4A006C, 0xD608309796947758, 0x397340ABDCD96956, 0x424B5DAE0CB63784, 0x000003726280F304 }; // XRB1
// Eve's generator values {XPC0 + XPC1*i, XQC0, XRC0 + XRC1*i} in GF(p747^2), expressed in normal representation
const uint64_t C_gen[5 * NWORDS64_FIELD]         = { 0xAF69BDDEC9296070, 0x8AC431344B2286BD, 0x3CFA47D203F07AFE, 0x162A8F46E4813F07, 0xAD4DDD2B67753675, 0x0E2EC4FDA5C93F08,
													 0xA676A39D0B8F01A0, 0xF5ED1D43A66A18AE, 0xA435E81C4D0EB5BB, 0x6CA414465FE77EB5, 0xAA8EB4A039EC4B7D, 0x000004566C7095BA,   // XPC0
													 0xB38034500C6DA1D2, 0x8F6EC8D9A1F35F28, 0xF8929FCCF0E08F28, 0xE26173136E9C4823, 0x40FCFEF0D82BE6AD, 0xD250DB7DCD87DA8A,
													 0x5D8128D2003719D1, 0xDD15896DE5C7F0EE, 0xE3A5A817AABA93DB, 0xB9A7EBF341C79B6E, 0x36057976E121CFDC, 0x0000010EABEEFEA0,   // XPC1
													 0xC9D02733A27AB49A, 0xB469BD77E0168E33, 0x05F8C5398CDFFBC7, 0x3E4A8125875936D8, 0x992DD94A7FF49581, 0x43A3E31079E1E5B6,
													 0x3E2A56DB507C88DB, 0xD066713B82EE0EA2, 0x0297C0C5A50BCB01, 0xD56B23D0DBB84C26, 0xE4E05108CB45392C, 0x000001C41F266159,   // XQC0
													 0x4A7E2CAF8075DBAE, 0x7C8CE9CE3F662D39, 0x0E5F171AAAD4D525, 0xD49B0EB806B01748, 0x6BD4262EE20D91E4, 0x8E0D5B740520C4D6,
													 0x1E04229F62707182, 0xF158168ED5A1579E, 0x69BAC9B55573B8F4, 0xB9FC03653052FCD9, 0x335A33155EB8B3ED, 0x000004830A950BFA,   // XRC0
													 0xB104DBDA0485994F, 0xC3AA33731C632A2E, 0x7BB8CEE8B3D9982B, 0xC430A10219BAF350, 0x67093EA63B360D7E, 0xF0FE015252925652,
													 0x5307546E0239541D, 0xA287B3C86C8687F2, 0xEDDF662A8E15DAB3, 0x49C23F9F35F33A30, 0xDBD16176640E8A0C, 0x0000014BE3B9B788 }; // XRC1

// Montgomery constant Montgomery_R2 = (2^768)^2 mod p747
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0xB72C5563CFD0070C, 0x3DCBDB82AE4B2EFB, 0x53B6DF3D0115B350, 0xED5F4AC6190451CD, 0xCF11EABCFB4DBBA5 ,0x4723FDABFDEA5C88, 
                                                     0x909485CA107DA103, 0xD233A15F550C0A86, 0x9F5B5C3A9349160E, 0xE1EC1C1C9606CD3C, 0x8E923055349AF253, 0x0000022019A271D0 };                                                    
// Value one in Montgomery representation 
const uint64_t Montgomery_one[NWORDS64_FIELD]    = { 0x000000000033EC27, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x1ED05A8A77BC4770, 0xCE23A20979D1B825,
                                                     0x7181DEF5FF042781, 0x69C3F4ABFB5A29DA, 0x890D3B44E106D57E, 0x3AE49E582C13F94E, 0xEA4A56D1578BCD2E, 0x000003A5F21C71B5 };
// x-coordinate of alpha which is a point with order 2 on the base curve = (i,0) -> XAlpha0 = 0 and XAlpha1 = 1
const uint64_t E0_alpha[NWORDS64_FIELD]          = { 0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
													 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, }; //xAlpha1
// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice] = { 
0, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 5, 6, 7, 8, 8, 9, 9, 9, 9,
9, 9, 9, 12, 11, 12, 12, 13, 14, 15, 16, 16, 16, 16, 16, 16, 17, 17, 18, 18, 17,
21, 17, 18, 21, 20, 21, 21, 21, 21, 21, 22, 25, 25, 25, 26, 27, 28, 28, 29, 30,
31, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 35, 36, 36, 33, 36, 35, 36, 36, 35,
36, 36, 37, 38, 38, 39, 40, 41, 42, 38, 39, 40, 41, 42, 40, 46, 42, 43, 46, 46,
46, 46, 48, 48, 48, 48, 49, 49, 48, 53, 54, 51, 52, 53, 54, 55, 56, 57, 58, 59,
59, 60, 62, 62, 63, 64, 64, 64 };

const unsigned int strat_Bob[MAX_Bob] = { 
0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10,
12, 12, 12, 12, 12, 12, 13, 14, 14, 15, 16, 16, 16, 16, 16, 17, 16, 16, 17, 19,
19, 20, 21, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 24, 24, 25, 27, 27, 28, 28,
29, 28, 29, 28, 28, 28, 30, 28, 28, 28, 29, 30, 33, 33, 33, 33, 34, 35, 37, 37,
37, 37, 38, 38, 37, 38, 38, 38, 38, 38, 39, 43, 38, 38, 38, 38, 43, 40, 41, 42,
43, 48, 45, 46, 47, 47, 48, 49, 49, 49, 50, 51, 50, 49, 49, 49, 49, 51, 49, 53,
50, 51, 50, 51, 51, 51, 52, 55, 55, 55, 56, 56, 56, 56, 56, 58, 58, 61, 61, 61,
63, 63, 63, 64, 65, 65, 65 };

const unsigned int strat_Eve[MAX_Eve] = { 
0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 3, 4, 4, 5, 5, 6, 5, 6, 6, 6, 7,8, 8, 9, 9, 9, 9,
9, 9, 9, 12, 10, 12, 12, 12, 12, 13, 12, 13, 13, 13, 14, 14, 14, 14, 18, 14, 18,
15, 17, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 20, 21,	22, 22, 22, 22, 23, 
23, 26, 23, 26, 23, 23, 26, 24, 26, 26, 27, 28, 27, 27, 28,	27, 28, 27, 28, 28,
28, 28, 29, 29, 31, 31, 31, 34, 34, 34, 34, 34, 34, 34, 34,	34, 34 };

// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions

#define fpcopy                  fpcopy747
#define fpzero                  fpzero747
#define fpadd                   fpadd747
#define fpsub                   fpsub747
#define fpneg                   fpneg747
#define fpdiv2                  fpdiv2_747
#define fpcorrection            fpcorrection747
#define fpmul_mont              fpmul747_mont
#define fpsqr_mont              fpsqr747_mont
#define fpinv_mont              fpinv747_mont
#define fpinv_chain_mont        fpinv747_chain_mont
#define fpinv_mont_bingcd       fpinv747_mont_bingcd
#define fp2copy                 fp2copy747
#define fp2zero                 fp2zero747
#define fp2add                  fp2add747
#define fp2sub                  fp2sub747
#define fp2neg                  fp2neg747
#define fp2div2                 fp2div2_747
#define fp2correction           fp2correction747
#define fp2mul_mont             fp2mul747_mont
#define fp2sqr_mont             fp2sqr747_mont
#define fp2inv_mont             fp2inv747_mont
#define fp2inv_mont_bingcd      fp2inv747_mont_bingcd
#define fpequal_non_constant_time  fpequal747_non_constant_time
#define mp_add_asm              mp_add747_asm
#define mp_addx2_asm            mp_add747x2_asm
#define mp_subx2_asm            mp_sub747x2_asm

#include "fpx.c"
#include "ec_isogeny.c"
#include "groupKey.c"
