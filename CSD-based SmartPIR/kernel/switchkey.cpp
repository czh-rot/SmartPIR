#include "top.h"

void intt_keyswitch_1_para1(
	UDTYPE target_in0[BRAMNUM][BRAMSIZE], UDTYPE target_out0[BRAMNUM][BRAMSIZE],
	UDTYPE idt0[N], UDTYPE sidt0[N], UDTYPE modulus0)
{
#pragma HLS INLINE off

	for (int j = 0; j < BRAMSIZE; j++)
	{
#pragma HLS PIPELINE
		for (int i = 0; i < BRAMNUM; i++)
		{
			target_out0[i][j] = target_in0[i][j];
		}
	}

	intt_4core<CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX>(
		target_out0, idt0, sidt0, modulus0);
}

void dyadic_mc1(UDTYPE buffer[BRAMNUM][BRAMSIZE],
				UDTYPEin *key_in0,
				UDTYPE buffer_out[2][BRAMNUM][BRAMSIZE],
				UDTYPE modulus, UDTYPE modulus_ratio0, UDTYPE modulus_ratio1)
{
#pragma HLS INLINE off

	UDTYPE2 buffer_temp[5];
#pragma HLS ARRAY_PARTITION variable = buffer_temp complete dim = 1

	UDTYPE2 buffer_temp_sum;
	UDTYPE local_wide_product[2];
#pragma HLS ARRAY_PARTITION variable = local_wide_product complete dim = 1
	UDTYPEin key_temp0, key_temp1;
	UDTYPE key0, key1, key2, key3, key4;
	UDTYPE out_data;

	IDXTYPE i, j, k;
	for (IDXTYPE idx = 0; idx < NN; idx++)
	{
#pragma HLS PIPELINE II = 1
#pragma HLS DEPENDENCE variable = buffer_out inter false

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);
		k = idx >> L_N;

		key_temp0 = key_in0[idx];

		key0 = key_temp0((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		buffer_temp[0] = buffer[i][j] * key0;

		out_data = buffer_out[k][i][j];

		buffer_temp_sum = (buffer_temp[0] + out_data);
		local_wide_product[0] = (buffer_temp_sum)(BITWIDTH - 1, 0);
		local_wide_product[1] = (buffer_temp_sum)(2 * BITWIDTH - 1, BITWIDTH);

		buffer_out[k][i][j] = barrett_reduce_128(local_wide_product[0], local_wide_product[1], modulus, modulus_ratio0, modulus_ratio1);
	}
}

void multiply_sub_add_cc_mc1(
	UDTYPEin *acc_sum_in, UDTYPEin *encrypted, UDTYPEin *out,
	UDTYPE poly_last[2][BRAMNUM][BRAMSIZE], UDTYPE poly[2][BRAMNUM][BRAMSIZE],
	UDTYPE modulus, UDTYPE ratio0, UDTYPE ratio1,
	UDTYPE modswitch_factors)
{

#pragma HLS INLINE off

	UDTYPEin encrypted_temp, acc_sum_temp, out_temp;
	UDTYPE encrypted_k0, encrypted_k1;
	UDTYPE acc_k0, acc_k1;

	UDTYPE encrypted_save_k0, encrypted_save_k1;
	UDTYPE temp_k0, temp_k1;

	IDXTYPE i, j;
	for (IDXTYPE idx = 0; idx < N; idx++)
	{
#pragma HLS PIPELINE II = 1
		//#pragma HLS DEPENDENCE variable = buffer_out inter false

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);

		// data read
		encrypted_temp = encrypted[idx];
		acc_sum_temp = acc_sum_in[idx];

		encrypted_k0 = encrypted_temp((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		encrypted_k1 = encrypted_temp((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);

		acc_k0 = acc_sum_temp((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		acc_k1 = acc_sum_temp((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);

		// k0
		temp_k0 = Hsub_single(poly[0][i][j], poly_last[0][i][j], modulus);
		temp_k0 = HMult_single_return(temp_k0, modswitch_factors, modulus, ratio0, ratio1);
		encrypted_save_k0 = Hadd_single(temp_k0, encrypted_k0, modulus);
		acc_k0 = Hadd_single(acc_k0, encrypted_save_k0, modulus);

		// k1
		temp_k1 = Hsub_single(poly[1][i][j], poly_last[1][i][j], modulus);
		temp_k1 = HMult_single_return(temp_k1, modswitch_factors, modulus, ratio0, ratio1);
		encrypted_save_k1 = Hadd_single(temp_k1, encrypted_k1, modulus);
		acc_k1 = Hadd_single(acc_k1, encrypted_save_k1, modulus);

		// data out
		out_temp((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH) = acc_k0;
		out_temp((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH) = acc_k1;

		out[idx] = out_temp;
	}
}

void data_out_mc1(UDTYPE in[2][BRAMNUM][BRAMSIZE], UDTYPEin *data_out0)
{

	UDTYPEin buffer_temp;

	IDXTYPE i, j;
	for (IDXTYPE idx = 0; idx < N; idx++)
	{
#pragma HLS PIPELINE II = 1
#pragma HLS DEPENDENCE variable = in inter false

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);

		buffer_temp((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH) = in[0][i][j];
		buffer_temp((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH) = in[1][i][j];

		data_out0[idx] = buffer_temp;
	}
}

void data_out_mc3(UDTYPE in[3][2][BRAMNUM][BRAMSIZE], UDTYPEin *data_out)
{

	UDTYPEin data_temp;
	IDXTYPE i, j, k;
	for (IDXTYPE idx = 0; idx < NN; idx++)
	{
#pragma HLS PIPELINE

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);
		k = idx >> L_N;

		(data_temp)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH) = in[0][k][i][j];
		(data_temp)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH) = in[1][k][i][j];
		(data_temp)((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH) = in[2][k][i][j];

		data_out[idx] = data_temp;
	}
}


void data_in_mc4(UDTYPEin *in,
				 UDTYPE operand[4][2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off
	UDTYPEin in_temp0, in_temp1;
	UDTYPEin op_tempk0_m0, op_tempk0_m1, op_tempk0_m2, op_tempk0_m3;
	UDTYPEin op_tempk1_m0, op_tempk1_m1, op_tempk1_m2;

	IDXTYPE i, j, k;
copy_data_in:
	for (IDXTYPE idx = 0; idx < NN; idx++)
	{
#pragma HLS PIPELINE
		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);
		k = idx >> L_N;

		in_temp0 = in[idx];

		op_tempk0_m0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		op_tempk0_m1 = (in_temp0)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		op_tempk0_m2 = (in_temp0)((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);
		op_tempk0_m3 = (in_temp0)((3 + 1) * BITWIDTH - 1, 3 * BITWIDTH);

		operand[0][k][i][j] = op_tempk0_m0;
		operand[1][k][i][j] = op_tempk0_m1;
		operand[2][k][i][j] = op_tempk0_m2;
		operand[3][k][i][j] = op_tempk0_m3;
	}
}

void pdata_in_mc4(UDTYPEin *in,
				 UDTYPE operand[4][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off
	UDTYPEin in_temp0, in_temp1;
	UDTYPEin op_tempk0_m0, op_tempk0_m1, op_tempk0_m2, op_tempk0_m3;
	UDTYPEin op_tempk1_m0, op_tempk1_m1, op_tempk1_m2;

	IDXTYPE i, j, k;
copy_data_in:
	for (IDXTYPE idx = 0; idx < N; idx++)
	{
#pragma HLS PIPELINE
		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);
//		k = idx >> L_N;

		in_temp0 = in[idx];

		op_tempk0_m0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		op_tempk0_m1 = (in_temp0)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		op_tempk0_m2 = (in_temp0)((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);
		op_tempk0_m3 = (in_temp0)((3 + 1) * BITWIDTH - 1, 3 * BITWIDTH);

		operand[0][i][j] = op_tempk0_m0;
		operand[1][i][j] = op_tempk0_m1;
		operand[2][i][j] = op_tempk0_m2;
		operand[3][i][j] = op_tempk0_m3;
	}
}

// in_k0��128bit�� ���������պ�����ģ��16384 cycles
void data_in_mc3(UDTYPEin *in,
				 UDTYPE operand[3][2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off
	UDTYPEin in_temp0, in_temp1;
	UDTYPEin op_tempk0_m0, op_tempk0_m1, op_tempk0_m2;
	UDTYPEin op_tempk1_m0, op_tempk1_m1, op_tempk1_m2;

	IDXTYPE i, j, k;
copy_data_in:
	for (IDXTYPE idx = 0; idx < NN; idx++)
	{
#pragma HLS PIPELINE
		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);
		k = idx >> L_N;

		in_temp0 = in[idx];

		op_tempk0_m0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		op_tempk0_m1 = (in_temp0)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		op_tempk0_m2 = (in_temp0)((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);

		operand[0][k][i][j] = op_tempk0_m0;
		operand[1][k][i][j] = op_tempk0_m1;
		operand[2][k][i][j] = op_tempk0_m2;
	}
}

void data_in_mc1(UDTYPEin *in,
				 UDTYPE operand[2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off
	UDTYPEin in_temp0, in_temp1;
	UDTYPEin op_tempk0_m0, op_tempk0_m1, op_tempk0_m2;
	UDTYPEin op_tempk1_m0, op_tempk1_m1, op_tempk1_m2;

	IDXTYPE i, j, k;
copy_data_in:
	for (IDXTYPE idx = 0; idx < NN; idx++)
	{
#pragma HLS PIPELINE
		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);
		k = idx >> L_N;

		in_temp0 = in[idx];

		op_tempk0_m0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);

		operand[k][i][j] = op_tempk0_m0;
	}
}
