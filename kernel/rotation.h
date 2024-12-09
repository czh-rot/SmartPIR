#include "define.h"

void input_galois_ntt_mc3(UDTYPEin *in_k0,
	UDTYPE operand[3][2][BRAMNUM][BRAMSIZE]);

void input_galois_ntt_mc1(UDTYPEin *in_k0,
		UDTYPE operand[2][BRAMNUM][BRAMSIZE]);

size_t reverse_bits(size_t operand, int bit_count);
void apply_galois_ntt(UDTYPE input[BRAMNUM][BRAMSIZE],
					  UDTYPE result[BRAMNUM][BRAMSIZE], size_t galois_elt);

void apply_galois_ntt_mc3(UDTYPE input0[BRAMNUM][BRAMSIZE], UDTYPE input1[BRAMNUM][BRAMSIZE],
		UDTYPE input2[BRAMNUM][BRAMSIZE], UDTYPE result0[BRAMNUM][BRAMSIZE],
		UDTYPE result1[BRAMNUM][BRAMSIZE], UDTYPE result2[BRAMNUM][BRAMSIZE],
					  size_t galois_elt);
void galois_inplace_mc3(UDTYPE operand[3][2][BRAMNUM][BRAMSIZE],
		UDTYPE target[3][BRAMNUM][BRAMSIZE], size_t galois_elt);

void galois_inplace_mc1(UDTYPE operand[1][2][BRAMNUM][BRAMSIZE],
		UDTYPE target[1][BRAMNUM][BRAMSIZE], size_t galois_elt);

void galois_mc1(
		UDTYPE in[2][BRAMNUM][BRAMSIZE],
		UDTYPE out[2][BRAMNUM][BRAMSIZE], UDTYPE target[BRAMNUM][BRAMSIZE],
		size_t galois_elt);

/*void apply_galois_inplace(
	UDTYPE in0[BRAMNUM][BRAMSIZE], UDTYPE in1[BRAMNUM][BRAMSIZE],
	UDTYPE out[2][BRAMNUM][BRAMSIZE], UDTYPE out2[BRAMNUM][BRAMSIZE],
	size_t galois_elt);*/



/*
template<unsigned mod_count>
void apply_galois_inplace_new(
	UDTYPE in[mod_count][2][BRAMNUM][BRAMSIZE],
	UDTYPE out[mod_count][2][BRAMNUM][BRAMSIZE], UDTYPE out2[mod_count][BRAMNUM][BRAMSIZE],
	size_t galois_elt)
{

	galois_loop:
	for (int i = 0; i < mod_count; i++){
#pragma HLS UNROLL
		apply_galois_inplace(in[i][0], in[i][1], out[i], out2[i], galois_elt);
	}
}
*/


