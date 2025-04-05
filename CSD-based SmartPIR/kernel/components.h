#pragma once
#include "define.h"

//--------------22-03-09
UDTYPE add_uint64(
    UDTYPE operand1, UDTYPE operand2, UDTYPE* result);

UDTYPE add_uint64_generic(
    UDTYPE operand1, UDTYPE operand2, uint8_t carry,
    UDTYPE* result);

UDTYPE Hadd_single(UDTYPE operand1,UDTYPE operand2,UDTYPE modulus_value);

uint8_t sub_uint64(
    UDTYPE operand1, UDTYPE operand2, UDTYPE* result);

uint8_t sub_uint64_generic(
    UDTYPE operand1, UDTYPE operand2,
    uint8_t borrow, UDTYPE* result);

UDTYPE sub_uint_uint_mod(
    UDTYPE operand1, UDTYPE operand2,
    UDTYPE modulus);

UDTYPE Hsub_single(UDTYPE operand1,UDTYPE operand2,UDTYPE modulus);

UDTYPE barrett_reduce_63(
    UDTYPE input, UDTYPE modulus, UDTYPE modulus_ratio1) ;

UDTYPE barrett_reduce_128(
    UDTYPE input0,UDTYPE input1, UDTYPE modulus, UDTYPE modulus_ratio0, UDTYPE modulus_ratio1);

UDTYPE HMult_single_return(UDTYPE operand1,UDTYPE operand2,
    UDTYPE modulus,UDTYPE modulus_ratio_0, UDTYPE modulus_ratio_1);

//--------------22-03-09

/*void Hadd_N (UDTYPE operand1[BRAMNUM][BRAMSIZE], UDTYPE operand2[BRAMNUM][BRAMSIZE], UDTYPE modulus);
UDTYPE Hadd_single(UDTYPE operand1,UDTYPE operand2,UDTYPE modulus_value);

UDTYPE HMult_single(UDTYPE operand1,UDTYPE operand2,
    UDTYPE modulus,UDTYPE modulus_ratio[2]);


// operand1 + operand2 => operand1
template<unsigned modcount>
void Hadd(UDTYPE operand1[modcount][2][BRAMNUM][BRAMSIZE],
	    UDTYPE operand2[modcount][2][BRAMNUM][BRAMSIZE],
	    UDTYPE modulus[INITMODCOUNT])
{
	for(int i = 0; i < modcount; i++){
#pragma HLS UNROLL
		for(int j = 0; j < 2; j++){
			Hadd_N(operand1[i][j], operand2[i][j], modulus[i]);
		}
	}
}


void multiply_uint64_hw64(
    UDTYPE operand1, UDTYPE operand2, UDTYPE* hw64);
void multiply_uint64(
    UDTYPE operand1, UDTYPE operand2, UDTYPE result128[2]);
void add_poly_poly_coeffmod(const UDTYPE* operand1,
    const UDTYPE* operand2,
    const UDTYPE modulus, UDTYPE* result);
uint8_t sub_uint64(
    UDTYPE operand1, UDTYPE operand2, UDTYPE* result);
void sub_poly_poly_coeffmod(const UDTYPE* operand1,
    const UDTYPE* operand2, size_t coeff_count,
    const UDTYPE modulus, UDTYPE* result);
uint8_t sub_uint64_generic(
    UDTYPE operand1, UDTYPE operand2,
    uint8_t borrow, UDTYPE* result);
UDTYPE sub_uint_uint_mod(
    UDTYPE operand1, UDTYPE operand2,
    const UDTYPE modulus);
UDTYPE barrett_reduce_63(
    UDTYPE input, const UDTYPE modulus, const UDTYPE modulus_ratio1);   //这个input其实可以63-bit�?


void barrett_reduce_N(UDTYPE operand1[BRAMNUM][BRAMSIZE],
	    UDTYPE modulus,UDTYPE modulus_ratio1);


void multiply_poly_scalar_coeffmod(const UDTYPE* poly,
    size_t coeff_count, UDTYPE scalar, const UDTYPE modulus, const UDTYPE modulus_ratio0, const UDTYPE modulus_ratio1,
    UDTYPE* result);
UDTYPE add_uint64(
    UDTYPE operand1, UDTYPE operand2, UDTYPE* result);
UDTYPE add_uint64_generic(
    UDTYPE operand1, UDTYPE operand2, uint8_t carry,
    UDTYPE* result);
void sub_poly_poly_coeffmod(const UDTYPE* operand1,
    const UDTYPE* operand2, size_t coeff_count,
    const UDTYPE modulus, UDTYPE* result);
UDTYPE barrett_reduce_128(
    UDTYPE* input, const UDTYPE modulus, UDTYPE modulus_ratio0, UDTYPE modulus_ratio1);
UDTYPE barrett_reduce_128(
    UDTYPE input0,UDTYPE input1, const UDTYPE modulus, UDTYPE modulus_ratio0, UDTYPE modulus_ratio1);

void multiply_coeffmod(UDTYPE operand1,UDTYPE operand2,
    const UDTYPE modulus, const UDTYPE modulus_ratio0, const UDTYPE modulus_ratio1,
    UDTYPE* result);

UDTYPE multiply_coeffmod(UDTYPE operand1,UDTYPE operand2,
    const UDTYPE modulus, const UDTYPE modulus_ratio0, const UDTYPE modulus_ratio1);

void add_8_1024(const UDTYPE operand1[BRAMNUM][BRAMSIZE],
    const UDTYPE operand2[BRAMNUM][BRAMSIZE],
    const UDTYPE modulus, UDTYPE result[BRAMNUM][BRAMSIZE]);

void sub_8_1024(const UDTYPE operand1[BRAMNUM][BRAMSIZE],
    const UDTYPE operand2[BRAMNUM][BRAMSIZE],
    const UDTYPE modulus, UDTYPE result[BRAMNUM][BRAMSIZE]);

void multiply_scalar_8_1024(const UDTYPE poly[BRAMNUM][BRAMSIZE],
    size_t coeff_count, UDTYPE scalar, const UDTYPE modulus, const UDTYPE modulus_ratio0, const UDTYPE modulus_ratio1,
    UDTYPE result[BRAMNUM][BRAMSIZE]);

UDTYPE HMult_single_return(UDTYPE operand1,UDTYPE operand2,
    UDTYPE modulus,UDTYPE modulus_ratio_0, UDTYPE modulus_ratio_1);

void HMult_single(UDTYPE operand1,UDTYPE operand2, UDTYPE &result,
    UDTYPE modulus,UDTYPE modulus_ratio[2]);

UDTYPE Hsub_single(UDTYPE operand1,UDTYPE operand2,UDTYPE modulus);*/
