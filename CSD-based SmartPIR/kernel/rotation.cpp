#include "top.h"

size_t reverse_bits(size_t operand, int bit_count)
{
	uint32_t res = 0;
	for (int i = 0; i < bit_count; ++i)
	{
		res |= ((operand >> i) & 1) << (bit_count - 1 - i);
	}
	return res;
}

void apply_galois_ntt(UDTYPE input[BRAMNUM][BRAMSIZE],
					  UDTYPE result[BRAMNUM][BRAMSIZE], size_t galois_elt)
{
	size_t coeff_count = size_t(1) << STAGENUM;
	size_t m_minus_one = 2 * coeff_count - 1;

	apply_i:
	for (size_t i = 0; i < BRAMNUM; i++)
	{
		apply_j:
		for (size_t j = 0; j < BRAMSIZE; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = input inter false
			size_t bits = i * BRAMSIZE + j;
			size_t reversed = reverse_bits(bits, STAGENUM);
			size_t index_raw = galois_elt * (2 * reversed + 1);
			index_raw &= m_minus_one;
			size_t index = reverse_bits((index_raw - 1) >> 1, STAGENUM);
			result[i][j] = input[index % BRAMNUM][index / BRAMNUM];
			input[index % BRAMNUM][index / BRAMNUM] = 0;
		}
	}
}



// in_k0是128bit， 三个数（刚好三个模）16384 cycles
void input_galois_ntt_mc3(UDTYPEin *in_k0,
		UDTYPE operand[3][2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off
	UDTYPEin in_temp0, in_temp1;
	UDTYPEin op_tempk0_m0, op_tempk0_m1, op_tempk0_m2;
	UDTYPEin op_tempk1_m0, op_tempk1_m1, op_tempk1_m2;

	/*copy_data_in:
	for (size_t k = 0; k < 2; k++){
		for (size_t j = 0; j < BRAMSIZE; j++){
			for (size_t i = 0; i < BRAMNUM; i++){
	#pragma HLS PIPELINE
				in_temp0 = *in_k0;
				in_k0++;

				op_tempk0_m0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
				op_tempk0_m1 = (in_temp0)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
				op_tempk0_m2 = (in_temp0)((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);

				operand[0][k][i][j] = op_tempk0_m0;
				operand[1][k][i][j] = op_tempk0_m1;
				operand[2][k][i][j] = op_tempk0_m2;
			}

		}
	}*/

	IDXTYPE i, j, k;
	copy_data_in:
	for (IDXTYPE idx = 0; idx < NN; idx++){
#pragma HLS PIPELINE
		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);
		k = idx >> L_N;

		in_temp0 = in_k0[idx];

		op_tempk0_m0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		op_tempk0_m1 = (in_temp0)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		op_tempk0_m2 = (in_temp0)((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);

		operand[0][k][i][j] = op_tempk0_m0;
		operand[1][k][i][j] = op_tempk0_m1;
		operand[2][k][i][j] = op_tempk0_m2;

	}
}


void input_galois_ntt_mc1(UDTYPEin *in_k0,
		UDTYPE operand[2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off
	UDTYPEin in_temp0, in_temp1;
	UDTYPEin op_tempk0_m0, op_tempk0_m1, op_tempk0_m2;
	UDTYPEin op_tempk1_m0, op_tempk1_m1, op_tempk1_m2;

	copy_data_in:
	for (size_t k = 0; k < 2; k++){
		for (size_t j = 0; j < BRAMSIZE; j++){
			for (size_t i = 0; i < BRAMNUM; i++){
#pragma HLS PIPELINE
				in_temp0 = *in_k0;
				in_k0++;

				op_tempk0_m0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
				operand[k][i][j] = op_tempk0_m0;

			}

		}
	}
}

void apply_galois_ntt_mc3(UDTYPE input0[BRAMNUM][BRAMSIZE], UDTYPE input1[BRAMNUM][BRAMSIZE],
		UDTYPE input2[BRAMNUM][BRAMSIZE], UDTYPE result0[BRAMNUM][BRAMSIZE],
		UDTYPE result1[BRAMNUM][BRAMSIZE], UDTYPE result2[BRAMNUM][BRAMSIZE],
					  size_t galois_elt)
{
	size_t coeff_count = size_t(1) << STAGENUM;
	size_t m_minus_one = 2 * coeff_count - 1;

	for (size_t i = 0; i < BRAMNUM; i++)
	{
		for (size_t j = 0; j < BRAMSIZE; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = input0 inter false
#pragma HLS DEPENDENCE variable = input1 inter false
#pragma HLS DEPENDENCE variable = input2 inter false

			size_t bits = i * BRAMSIZE + j;
			size_t reversed = reverse_bits(bits, STAGENUM);
			size_t index_raw = galois_elt * (2 * reversed + 1);
			index_raw &= m_minus_one;
			size_t index = reverse_bits((index_raw - 1) >> 1, STAGENUM);

			result0[i][j] = input0[index % BRAMNUM][index / BRAMNUM];
			input0[index % BRAMNUM][index / BRAMNUM] = 0;

			result1[i][j] = input1[index % BRAMNUM][index / BRAMNUM];
			input1[index % BRAMNUM][index / BRAMNUM] = 0;

			result2[i][j] = input2[index % BRAMNUM][index / BRAMNUM];
			input2[index % BRAMNUM][index / BRAMNUM] = 0;
		}
	}
}



void galois_inplace_mc3(UDTYPE operand[3][2][BRAMNUM][BRAMSIZE],
		UDTYPE target[3][BRAMNUM][BRAMSIZE], size_t galois_elt){
#pragma HLS INLINE off

	apply_galois_ntt_mc3(operand[0][0], operand[1][0], operand[2][0],
			 target[0], target[1],  target[2], galois_elt);

	apply_galois_ntt_mc3(operand[0][1], operand[1][1], operand[2][1],
			operand[0][0], operand[1][0], operand[2][0],  galois_elt);


	UDTYPE temp1;

	re_index:
	for (size_t j = 0; j < BRAMSIZE; j++){
#pragma HLS PIPELINE
		for (size_t i = 0; i < BRAMNUM; i++){
			for (size_t m = 0; m < 3; m++){
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = target inter false

				temp1 = operand[m][0][i][j];
				operand[m][0][i][j] = target[m][i][j];
				target[m][i][j] = temp1;

			}
		}
	}
}

void galois_inplace_mc1(UDTYPE operand[1][2][BRAMNUM][BRAMSIZE],
		UDTYPE target[1][BRAMNUM][BRAMSIZE], size_t galois_elt){
#pragma HLS INLINE off

	apply_galois_ntt(operand[0][0], target[0], galois_elt);

	apply_galois_ntt(operand[0][1], operand[0][0], galois_elt);


	UDTYPE temp1;

	re_index:
	for (size_t j = 0; j < BRAMSIZE; j++){
#pragma HLS PIPELINE
		for (size_t i = 0; i < BRAMNUM; i++){
			for (size_t m = 0; m < 1; m++){
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = target inter false

				temp1 = operand[m][0][i][j];
				operand[m][0][i][j] = target[m][i][j];
				target[m][i][j] = temp1;

			}
		}
	}
}

void galois_mc1(
		UDTYPE in[2][BRAMNUM][BRAMSIZE],
		UDTYPE out[2][BRAMNUM][BRAMSIZE], UDTYPE target[BRAMNUM][BRAMSIZE],
		size_t galois_elt)
{
#pragma HLS INLINE off

	apply_galois_ntt(in[0], out[0], galois_elt);
	apply_galois_ntt(in[1], target, galois_elt);

}

/*void apply_galois_inplace_mc3(UDTYPE operand[3][2][BRAMNUM][BRAMSIZE],
		UDTYPE target[3][BRAMNUM][BRAMSIZE], size_t galois_elt){
#pragma HLS INLINE off


	galois_loop:
	for (int i = 0; i < 3; i++){
#pragma HLS UNROLL
		apply_galois_ntt_double(operand[i], target[i], galois_elt);
		//apply_galois_ntt(operand[i][0], galois_elt, target[i]);
		//apply_galois_ntt(operand[i][1], galois_elt, operand[i][0]);
	}


	UDTYPE temp1;

	re_index:
	for (size_t j = 0; j < BRAMSIZE; j++){
#pragma HLS PIPELINE
		for (size_t i = 0; i < BRAMNUM; i++){
			for (size_t m = 0; m < 3; m++){
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = target inter false

				temp1 = operand[m][0][i][j];
				operand[m][0][i][j] = target[m][i][j];
				target[m][i][j] = temp1;

			}
		}
	}
}*/
