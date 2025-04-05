#include "define.h"
#include "rotation.h"

void data_out_mc1(UDTYPE in[2][BRAMNUM][BRAMSIZE], UDTYPEin *data_out0);
void data_out_mc3(UDTYPE in[3][2][BRAMNUM][BRAMSIZE], UDTYPEin *data_out);
void data_in_mc3(UDTYPEin *in, UDTYPE operand[3][2][BRAMNUM][BRAMSIZE]);
void data_in_mc4(UDTYPEin *in, UDTYPE operand[4][2][BRAMNUM][BRAMSIZE]);
void data_in_mc1(UDTYPEin *in, UDTYPE operand[2][BRAMNUM][BRAMSIZE]);
void pdata_in_mc4(UDTYPEin *in,UDTYPE operand[4][BRAMNUM][BRAMSIZE]);

void multiply_sub_add_cc_mc1(
	UDTYPEin *acc_sum_in, UDTYPEin *encrypted, UDTYPEin *out,
	UDTYPE poly_last[2][BRAMNUM][BRAMSIZE], UDTYPE poly[2][BRAMNUM][BRAMSIZE],
	UDTYPE modulus, UDTYPE ratio0, UDTYPE ratio1,
	UDTYPE modswitch_factors);

void intt_keyswitch_1_para1(
	UDTYPE target_in0[BRAMNUM][BRAMSIZE], UDTYPE target_out0[BRAMNUM][BRAMSIZE],
	UDTYPE idt0[N], UDTYPE sidt0[N], UDTYPE modulus0);

void dyadic_mc1(UDTYPE buffer[BRAMNUM][BRAMSIZE],
				UDTYPEin *key_in0,
				UDTYPE buffer_out[2][BRAMNUM][BRAMSIZE],
				UDTYPE modulus, UDTYPE modulus_ratio0, UDTYPE modulus_ratio1);

template <unsigned mod_count>
void galois(
	UDTYPE in[mod_count][2][BRAMNUM][BRAMSIZE],
	UDTYPE out[mod_count][2][BRAMNUM][BRAMSIZE], UDTYPE target[mod_count][BRAMNUM][BRAMSIZE],
	size_t galois_elt)
{
#pragma HLS INLINE off

galois_loop:
	for (int i = 0; i < mod_count; i++)
	{
#pragma HLS UNROLL
		apply_galois_ntt(in[i][0], out[i][0], galois_elt);
		apply_galois_ntt(in[i][1], target[i], galois_elt);
	}
}

// 这个代码之后可能会需要参�??-----------------------------------------------------------------------------
/*template <unsigned portnum, unsigned modcount, unsigned bramnum, unsigned bramsize>
void Dyadic_keyswitch_manymod(UDTYPE buffer[modcount][bramnum][bramsize],
							  UDTYPEin *key_in0, UDTYPEin *key_in1, UDTYPEin *key_in2, UDTYPEin *key_in3,
							  UDTYPE buffer_out[2][bramnum][bramsize],
							  UDTYPE modulus, UDTYPE modulus_ratio0, UDTYPE modulus_ratio1)
{
#pragma HLS INLINE off

	UDTYPE2 buffer_temp[modcount];
#pragma HLS ARRAY_PARTITION variable = buffer_temp complete dim = 1
	UDTYPE key[modcount];
#pragma HLS ARRAY_PARTITION variable = key complete dim = 1

	UDTYPE2 buffer_temp_sum;

	UDTYPE local_wide_product[2];
#pragma HLS ARRAY_PARTITION variable = local_wide_product complete dim = 1

	UDTYPEin key_temp;

	int index = 0;

	for (int k = 0; k < 2; k++)
	{
		for (int j = 0; j < bramsize; j++)
		{
			for (int n = 0; n < bramnum; n += portnum)
			{
#pragma HLS PIPELINE II = 1

			Dyadic_2_label0:
				for (int i = 0; i < portnum; i++)
				{ // PORTNUM可以是设�??1�??2�??3�??4�??

					if (i == 0)
					{
						key_temp = (key_in0[index]);
					}
					else if (i == 1)
					{
						key_temp = (key_in1[index]);
					}
					else if (i == 2)
					{
						key_temp = (key_in2[index]);
					}
					else if (i == 3)
					{
						key_temp = (key_in3[index]);
					}

					buffer_temp_sum = 0;

					for (int m = 0; m < modcount; m++)
					{
						key[m] = key_temp((m + 1) * BITWIDTH - 1, m * BITWIDTH);
						buffer_temp[m] = buffer[m][n + i][j] * key[m];
						buffer_temp_sum += buffer_temp[m];
					}

					local_wide_product[0] = (buffer_temp_sum)(BITWIDTH - 1, 0);
					local_wide_product[1] = (buffer_temp_sum)(2 * BITWIDTH - 1, BITWIDTH);
					buffer_out[k][n + i][j] = barrett_reduce_128(local_wide_product[0], local_wide_product[1], modulus, modulus_ratio0, modulus_ratio1);
				}
				index++;
			}
		}
	}
}*/
// -----------------------------------------------------------------------------

//当UDTYPE小于32bit时，modcount�??多为4
template <unsigned modcount>
void dyadic_unpara_poly(UDTYPE buffer[modcount][BRAMNUM][BRAMSIZE],
						UDTYPEin *key_in0,
						UDTYPE buffer_out[2][BRAMNUM][BRAMSIZE],
						UDTYPE modulus, UDTYPE modulus_ratio0, UDTYPE modulus_ratio1)
{
#pragma HLS INLINE off

	UDTYPE2 buffer_temp[modcount];
#pragma HLS ARRAY_PARTITION variable = buffer_temp complete dim = 1
	UDTYPE key[modcount];
#pragma HLS ARRAY_PARTITION variable = key complete dim = 1

	UDTYPE2 buffer_temp_sum;

	UDTYPE local_wide_product[2];
#pragma HLS ARRAY_PARTITION variable = local_wide_product complete dim = 1

	UDTYPEin key_temp;

	IDXTYPE i, j, k;
	for (IDXTYPE idx = 0; idx < NN; idx++)
	{
#pragma HLS PIPELINE II = 1

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);
		k = idx >> L_N;

		key_temp = key_in0[idx];
		buffer_temp_sum = 0;

		for (IDXTYPE m = 0; m < modcount; m++)
		{
			key[m] = key_temp((m + 1) * BITWIDTH - 1, m * BITWIDTH);
			buffer_temp[m] = buffer[m][i][j] * key[m];
			buffer_temp_sum += buffer_temp[m];
		}

		local_wide_product[0] = (buffer_temp_sum)(BITWIDTH - 1, 0);
		local_wide_product[1] = (buffer_temp_sum)(2 * BITWIDTH - 1, BITWIDTH);
		buffer_out[k][i][j] = barrett_reduce_128(local_wide_product[0], local_wide_product[1], modulus, modulus_ratio0, modulus_ratio1);
	}
}

template <unsigned modcount>
void dyadic_unpara_poly_new(UDTYPE buffer[modcount][BRAMNUM][BRAMSIZE],
							UDTYPEin *key_in0,
							UDTYPE buffer_out[2][BRAMNUM][BRAMSIZE],
							UDTYPE modulus, UDTYPE modulus_ratio0, UDTYPE modulus_ratio1,
							size_t index_j, UDTYPE target[BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	UDTYPE2 buffer_temp[modcount];
#pragma HLS ARRAY_PARTITION variable = buffer_temp complete dim = 1
	UDTYPE key[modcount];
#pragma HLS ARRAY_PARTITION variable = key complete dim = 1

	UDTYPE2 buffer_temp_sum;

	UDTYPE local_wide_product[2];
#pragma HLS ARRAY_PARTITION variable = local_wide_product complete dim = 1

	UDTYPEin key_temp;

	IDXTYPE i, j, k;
	for (IDXTYPE idx = 0; idx < NN; idx++)
	{
#pragma HLS PIPELINE II = 1

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);
		k = idx >> L_N;

		key_temp = key_in0[idx];
		buffer_temp_sum = 0;

		for (IDXTYPE m = 0; m < modcount; m++)
		{
			key[m] = key_temp((m + 1) * BITWIDTH - 1, m * BITWIDTH);

			if (m == index_j)
			{
				buffer_temp[m] = target[i][j] * key[m];
			}
			else
			{
				buffer_temp[m] = buffer[m][i][j] * key[m];
			}

			buffer_temp_sum += buffer_temp[m];
		}

		local_wide_product[0] = (buffer_temp_sum)(BITWIDTH - 1, 0);
		local_wide_product[1] = (buffer_temp_sum)(2 * BITWIDTH - 1, BITWIDTH);
		buffer_out[k][i][j] = barrett_reduce_128(local_wide_product[0], local_wide_product[1], modulus, modulus_ratio0, modulus_ratio1);
	}
}

template <unsigned corenum>
void multiply_sub(UDTYPE poly_last[2][BRAMNUM][BRAMSIZE], UDTYPE poly[2][BRAMNUM][BRAMSIZE],
				  UDTYPE encrypted[2][BRAMNUM][BRAMSIZE], UDTYPE modulus, UDTYPE ratio0, UDTYPE ratio1,
				  UDTYPE modswitch_factors)
{

#pragma HLS INLINE off

	for (size_t j = 0; j < BRAMSIZE; j++)
	{
		for (size_t m = 0; m < 2; m++)
		{
			for (size_t i = 0; i < BRAMNUM; i += corenum)
			{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = encrypted inter false

				for (size_t n = 0; n < corenum; n++)
				{

					size_t i_n;
					i_n = i + n;

					UDTYPE temp;
					temp = Hsub_single(poly[m][i_n][j], poly_last[m][i_n][j], modulus);
					temp = HMult_single_return(temp, modswitch_factors, modulus, ratio0, ratio1);
					encrypted[m][i_n][j] = Hadd_single(temp, encrypted[m][i_n][j], modulus);

					/*UDTYPE temp_result1;
					DTYPE borrow1 = sub_uint64(poly[m][i_n][j], poly_last[m][i_n][j], &temp_result1);
					temp1 = temp_result1 + (modulus & UDTYPE(-borrow1));

					temp1 = multiply_coeffmod(temp1, modswitch_factors,
											  modulus, modulus_ratio0, modulus_ratio1);

					UDTYPE sum = temp1 + encrypted[m][i_n][j];
					encrypted[m][i_n][j] = sum - (modulus & UDTYPE(-DTYPE(sum >= modulus)));*/
				}
			}
		}
	}
}

template <unsigned modcount>
void copy_intt(UDTYPE target_in[modcount][BRAMNUM][BRAMSIZE], UDTYPE target_out[modcount][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	for (int j = 0; j < BRAMSIZE; j++)
	{
#pragma HLS PIPELINE
		for (int m = 0; m < modcount; m++)
		{
			for (int i = 0; i < BRAMNUM; i++)
			{
				target_out[m][i][j] = target_in[m][i][j];
			}
		}
	}
}

// template <unsigned modcount, unsigned barcount>
// void ntt_keyswitch_1(UDTYPE target[modcount][BRAMNUM][BRAMSIZE], UDTYPE operand[modcount][BRAMNUM][BRAMSIZE],
// 					 UDTYPE rp[N], UDTYPE srp[N], UDTYPE modulus, UDTYPE modulus_ratio1)
// {
// #pragma HLS INLINE off

// 	UDTYPE operand_temp;

// ntt_j:
// 	for (int j = 0; j < BRAMSIZE; j++)
// 	{
// 	ntt_i:
// 		for (int i = 0; i < BRAMNUM; i += barcount)
// 		{
// #pragma HLS PIPELINE
// 			for (int c = 0; c < barcount; c++)
// 			{
// 				for (int m = 0; m < modcount; m++)
// 				{
// 					size_t i_c;
// 					i_c = i + c;

// 					operand_temp = target[m][i_c][j];
// 					operand[m][i_c][j] = barrett_reduce_63(operand_temp, modulus, modulus_ratio1);
// 				}
// 			}
// 		}
// 	}

// 	ntt_4core_mods_new<modcount, CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX, 0>(
// 		operand, rp, srp, modulus);
// }

template <unsigned modcount, unsigned barcount>
void ntt_keyswitch_1_new(UDTYPE target[modcount][BRAMNUM][BRAMSIZE], UDTYPE operand[modcount][BRAMNUM][BRAMSIZE],
						 UDTYPE rp[RPBRAMNUM][RPBRAMSIZE], UDTYPE srp[RPBRAMNUM][RPBRAMSIZE],
						 UDTYPE modulus_i[INITMODCOUNT], UDTYPE modulus, UDTYPE modulus_ratio1)
{
#pragma HLS INLINE off

	UDTYPE operand_temp;
	// size_t index = (j == modcount ? INITMODCOUNT - 1 : j);

ntt_j:
	for (int j = 0; j < BRAMSIZE; j++)
	{
	ntt_i:
		for (int i = 0; i < BRAMNUM; i += barcount)
		{
#pragma HLS PIPELINE
			for (int c = 0; c < barcount; c++)
			{
				for (int m = 0; m < modcount; m++)
				{

					size_t i_c;
					i_c = i + c;

					if (modulus_i[m] <= modulus)
					{
						operand[m][i_c][j] = target[m][i_c][j];
					}
					else
					{
						operand_temp = target[m][i_c][j];
						operand[m][i_c][j] = barrett_reduce_63(operand_temp, modulus, modulus_ratio1);
					}
				}
			}
		}
	}

	ntt_8core_mods_new<modcount, CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX, 0>(
		operand, rp, srp, modulus);
}

// template <unsigned para_i>
// void ntt_keyswitch_2(UDTYPE in[2][BRAMNUM][BRAMSIZE], UDTYPE out[2][BRAMNUM][BRAMSIZE],
// 					 UDTYPE rp[N], UDTYPE srp[N], UDTYPE modulus, UDTYPE modulus_ratio1, UDTYPE half)
// {
// #pragma HLS INLINE off

// 	UDTYPE operand_temp;
// 	UDTYPE half_mod = barrett_reduce_63(half, modulus, modulus_ratio1);

// 	for (size_t j = 0; j < BRAMSIZE; j++)
// 	{
// 		for (size_t i = 0; i < BRAMNUM; i += para_i)
// 		{
// #pragma HLS PIPELINE
// 			//#pragma HLS DEPENDENCE variable = in inter false
// 			for (size_t n = 0; n < para_i; n++)
// 			{
// 				for (size_t k = 0; k < 2; k++)
// 				{
// 					size_t i_n;
// 					i_n = i + n;

// 					operand_temp = barrett_reduce_63(in[k][i_n][j], modulus, modulus_ratio1);
// 					out[k][i_n][j] = sub_uint_uint_mod(operand_temp, half_mod, modulus);
// 				}
// 			}
// 		}
// 	}

// 	ntt_4core_mods<2, CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX, 1>(
// 		out, rp, srp, modulus);
// }

template <unsigned para_i>
void ntt_keyswitch_2_new(UDTYPE in[2][BRAMNUM][BRAMSIZE], UDTYPE out[2][BRAMNUM][BRAMSIZE],
						 UDTYPE rp[RPBRAMNUM][RPBRAMSIZE], UDTYPE srp[RPBRAMNUM][RPBRAMSIZE],
						 UDTYPE modulus, UDTYPE modulus_ratio1, UDTYPE half)
{
#pragma HLS INLINE off

	UDTYPE operand_temp;
	UDTYPE half_mod = barrett_reduce_63(half, modulus, modulus_ratio1);
	// cout << "half = " << half << endl;
	// cout << "half_mod = " << half_mod << endl;
	// cout << "q6-7-1 -- buffer_NTT1_p0 -- k0" << endl;

	for (size_t j = 0; j < BRAMSIZE; j++)
	{
		for (size_t i = 0; i < BRAMNUM; i += para_i)
		{
#pragma HLS PIPELINE
			//#pragma HLS DEPENDENCE variable = in inter false
			for (size_t n = 0; n < para_i; n++)
			{
				for (size_t k = 0; k < 2; k++)
				{
					size_t i_n;
					i_n = i + n;

					operand_temp = barrett_reduce_63(in[k][i_n][j], modulus, modulus_ratio1);
					// cout << "operand_temp = " << operand_temp << endl;
					out[k][i_n][j] = sub_uint_uint_mod(operand_temp, half_mod, modulus);
				}
			}
		}
	}

	// q6-7-2: j = 0, input of the NTT1
	// cout << "q6-7-2 -- buffer_NTT1_p0 -- k0" << endl;
	// cout << out[0][0][0] << endl;	 // 0
	// cout << out[0][1][0] << endl;	 // 1
	// cout << out[0][2][0] << endl;	 // 2
	// cout << out[0][3][0] << endl;	 // 3
	// cout << out[0][4][0] << endl;	 // 4
	// cout << out[0][5][0] << endl;	 // 5
	// cout << out[0][1][1] << endl;	 // 9
	// cout << out[0][0][512] << endl;	 // 4096
	// cout << out[0][0][900] << endl;	 // 7200
	// cout << out[0][1][900] << endl;	 // 7201
	// cout << out[0][7][1023] << endl; // 8191

	ntt_8core_mods_new<2, CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX, 1>(
		out, rp, srp, modulus);
}

// template <unsigned modcount>
// void intt_keyswitch_1(UDTYPE target_in[modcount][BRAMNUM][BRAMSIZE], UDTYPE target_out[modcount][BRAMNUM][BRAMSIZE],
// 					  UDTYPE idt[ROOTMODCOUNT][N], UDTYPE sidt[ROOTMODCOUNT][N], UDTYPE modulus[ROOTMODCOUNT])
// {
// #pragma HLS INLINE off

// 	copy_intt<modcount>(target_in, target_out);

// 	for (int i = 0; i < modcount; i++)
// 	{
// #pragma HLS UNROLL
// 		intt_4core_new<CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX, 1>(
// 			target_out[i], idt[i], sidt[i], modulus[i]);
// 	}
// }

template <unsigned modcount>
void intt_keyswitch_1_new(UDTYPE target_in[modcount][BRAMNUM][BRAMSIZE], UDTYPE target_out[modcount][BRAMNUM][BRAMSIZE],
						  UDTYPE idt[ROOTMODCOUNT][RPBRAMNUM][RPBRAMSIZE], UDTYPE sidt[ROOTMODCOUNT][RPBRAMNUM][RPBRAMSIZE],
						  UDTYPE modulus[ROOTMODCOUNT])
{
#pragma HLS INLINE off

	copy_intt<modcount>(target_in, target_out);

	for (int i = 0; i < modcount; i++)
	{
#pragma HLS UNROLL
		intt_8core_new<1, CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX, 1>(
			target_out[i], idt[i], sidt[i], modulus[i]);
	}
}

// template <unsigned para_i>
// void intt_keyswitch_2(UDTYPE operand[2][BRAMNUM][BRAMSIZE],
// 					  UDTYPE idt[N], UDTYPE sidt[N], UDTYPE modulus, UDTYPE modulus_ratio1)
// {
// #pragma HLS INLINE off

// 	intt_4core_mods_new<2, CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX, 0>(
// 		operand, idt, sidt, modulus);

// 	UDTYPE half = modulus >> 1;
// 	UDTYPE operand_temp;
// 	for (size_t j = 0; j < BRAMSIZE; j++)
// 	{
// 		for (size_t i = 0; i < BRAMNUM; i += para_i)
// 		{
// #pragma HLS PIPELINE
// #pragma HLS DEPENDENCE variable = operand inter false
// 			for (size_t n = 0; n < para_i; n++)
// 			{
// 				for (size_t k = 0; k < 2; k++)
// 				{
// 					size_t i_n;
// 					i_n = i + n;
// 					operand_temp = operand[k][i_n][j] + half;
// 					operand[k][i_n][j] = barrett_reduce_63(operand_temp, modulus, modulus_ratio1);
// 				}
// 			}
// 		}
// 	}
// }
template <unsigned para_i>
void intt_keyswitch_2_new(UDTYPE operand[2][BRAMNUM][BRAMSIZE],
						  UDTYPE idt[RPBRAMNUM][RPBRAMSIZE], UDTYPE sidt[RPBRAMNUM][RPBRAMSIZE],
						  UDTYPE modulus, UDTYPE modulus_ratio1)
{
#pragma HLS INLINE off

	intt_8core_mods_new<2, CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX, 0>(
		operand, idt, sidt, modulus);

	UDTYPE half = modulus >> 1;
	UDTYPE operand_temp;
	for (size_t j = 0; j < BRAMSIZE; j++)
	{
		for (size_t i = 0; i < BRAMNUM; i += para_i)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
			for (size_t n = 0; n < para_i; n++)
			{
				for (size_t k = 0; k < 2; k++)
				{
					size_t i_n;
					i_n = i + n;
					operand_temp = operand[k][i_n][j] + half;
					operand[k][i_n][j] = barrett_reduce_63(operand_temp, modulus, modulus_ratio1);
				}
			}
		}
	}
}

template <unsigned para_i>
void ntt_keyswitch_1_para1(UDTYPE target[BRAMNUM][BRAMSIZE], UDTYPE operand[BRAMNUM][BRAMSIZE],
						   UDTYPE rp[N], UDTYPE srp[N], UDTYPE modulus, UDTYPE modulus_ratio1)
{
#pragma HLS INLINE off

	UDTYPE operand_temp;

	for (int j = 0; j < BRAMSIZE; j++)
	{
		for (int i = 0; i < BRAMNUM; i += para_i)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = target inter false

			for (int c = 0; c < 2; c++)
			{
				size_t i_c;
				i_c = i + c;

				operand_temp = target[i_c][j];
				operand[i_c][j] = barrett_reduce_63(operand_temp, modulus, modulus_ratio1);
			}
		}
	}

	ntt_4core<CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX>(
		operand, rp, srp, modulus);
}
/*
template <unsigned modcount, unsigned ntt1_para, unsigned intt2_para, unsigned ntt2_para>
void switchkey_k5_para1_low(UDTYPE encrypted[modcount][2][BRAMNUM][BRAMSIZE],
							UDTYPE target[modcount][BRAMNUM][BRAMSIZE],
							UDTYPEin *keys_in0,
							UDTYPE rp[ROOTMODCOUNT][N],
							UDTYPE srp[ROOTMODCOUNT][N],
							UDTYPE idt[ROOTMODCOUNT][N],
							UDTYPE sidt[ROOTMODCOUNT][N],
							UDTYPE modulus[ROOTMODCOUNT],
							UDTYPE modulus_ratio[ROOTMODCOUNT][2],
							UDTYPE modswitch_factors[modcount],
							UDTYPE target_INTT[modcount][BRAMNUM][BRAMSIZE],
							UDTYPE buffer_NTT0_one_p0[BRAMNUM][BRAMSIZE],
							UDTYPE buffer_NTT1_p0[2][BRAMNUM][BRAMSIZE],
							UDTYPE buffer_dyadic_last[2][BRAMNUM][BRAMSIZE],
							UDTYPE buffer_dyadic_p0[2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	UDTYPE half = modulus[INITMODCOUNT - 1] >> 1;

	UDTYPE key_index = 0;

	intt_keyswitch_1_para1(target[0], target_INTT[0], idt[0], sidt[0], modulus[0]);
	intt_keyswitch_1_para1(target[1], target_INTT[1], idt[1], sidt[1], modulus[1]);
	intt_keyswitch_1_para1(target[2], target_INTT[2], idt[2], sidt[2], modulus[2]);
	intt_keyswitch_1_para1(target[3], target_INTT[3], idt[3], sidt[3], modulus[3]);
	intt_keyswitch_1_para1(target[4], target_INTT[4], idt[4], sidt[4], modulus[4]);

	int mod_idx = 5;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[0], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);

	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[1], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[2], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[5][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[3], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[4], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_2<intt2_para>(
		buffer_dyadic_last, idt[mod_idx], sidt[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	// 0

	mod_idx = 0;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[0], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][mod_idx], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[1], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[2], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[3], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[4], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

	// 1

	mod_idx = 1;
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[0], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[1], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[2], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[3], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[4], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

	// 2
	mod_idx = 2;
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[0], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[1], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[2], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[3], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[4], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

	// 3

	mod_idx = 3;
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[0], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[1], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[2], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[3], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[4], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

	// 4
	mod_idx = 4;
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[0], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[1], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[2], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[3], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT[4], buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);
}
*/

/*
template <unsigned modcount, unsigned ntt1_para, unsigned intt2_para, unsigned ntt2_para>
void switchkey_k5_para1(UDTYPE encrypted[modcount][2][BRAMNUM][BRAMSIZE],
						UDTYPE target[modcount][BRAMNUM][BRAMSIZE],
						UDTYPEin *keys_in0,
						UDTYPE rp[ROOTMODCOUNT][N],
						UDTYPE srp[ROOTMODCOUNT][N],
						UDTYPE idt[ROOTMODCOUNT][N],
						UDTYPE sidt[ROOTMODCOUNT][N],
						UDTYPE modulus[ROOTMODCOUNT],
						UDTYPE modulus_ratio[ROOTMODCOUNT][2],
						UDTYPE modswitch_factors[modcount],
						UDTYPE target_INTT_one_p0[BRAMNUM][BRAMSIZE],
						UDTYPE target_INTT_one_p1[BRAMNUM][BRAMSIZE],
						UDTYPE buffer_NTT0_one_p0[BRAMNUM][BRAMSIZE],
						UDTYPE buffer_NTT0_one_p1[BRAMNUM][BRAMSIZE],
						UDTYPE buffer_NTT1_p0[2][BRAMNUM][BRAMSIZE],
						UDTYPE buffer_NTT1_p1[2][BRAMNUM][BRAMSIZE],
						UDTYPE buffer_dyadic_last[2][BRAMNUM][BRAMSIZE],
						UDTYPE buffer_dyadic_p0[2][BRAMNUM][BRAMSIZE],
						UDTYPE buffer_dyadic_p1[2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE

	UDTYPE half = modulus[INITMODCOUNT - 1] >> 1;

	UDTYPE key_index = 0;

	int mod_idx = 5;

	intt_keyswitch_1_para1(target[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);

	key_index += 2 * N;

	intt_keyswitch_1_para1(target[1], target_INTT_one_p1, idt[1], sidt[1], modulus[1]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[2], target_INTT_one_p0, idt[2], sidt[2], modulus[2]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[5][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[3], target_INTT_one_p1, idt[3], sidt[3], modulus[3]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[4], target_INTT_one_p0, idt[4], sidt[4], modulus[4]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_2<intt2_para>(
		buffer_dyadic_last, idt[mod_idx], sidt[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	// 0

	mod_idx = 0;

	intt_keyswitch_1_para1(target[0], target_INTT_one_p1, idt[0], sidt[0], modulus[0]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][mod_idx], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[1], target_INTT_one_p0, idt[1], sidt[1], modulus[1]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[2], target_INTT_one_p1, idt[2], sidt[2], modulus[2]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[3], target_INTT_one_p0, idt[3], sidt[3], modulus[3]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[4], target_INTT_one_p1, idt[4], sidt[4], modulus[4]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

	// 1

	mod_idx = 1;

	intt_keyswitch_1_para1(target[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[1], target_INTT_one_p1, idt[1], sidt[1], modulus[1]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[2], target_INTT_one_p0, idt[2], sidt[2], modulus[2]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[3], target_INTT_one_p1, idt[3], sidt[3], modulus[3]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[4], target_INTT_one_p0, idt[4], sidt[4], modulus[4]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p1, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

	// 2
	mod_idx = 2;
	intt_keyswitch_1_para1(target[0], target_INTT_one_p1, idt[0], sidt[0], modulus[0]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[1], target_INTT_one_p0, idt[1], sidt[1], modulus[1]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[2], target_INTT_one_p1, idt[2], sidt[2], modulus[2]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[3], target_INTT_one_p0, idt[3], sidt[3], modulus[3]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[4], target_INTT_one_p1, idt[4], sidt[4], modulus[4]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

	// 3

	mod_idx = 3;
	intt_keyswitch_1_para1(target[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[1], target_INTT_one_p1, idt[1], sidt[1], modulus[1]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[2], target_INTT_one_p0, idt[2], sidt[2], modulus[2]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[3], target_INTT_one_p1, idt[3], sidt[3], modulus[3]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[4], target_INTT_one_p0, idt[4], sidt[4], modulus[4]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p1, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

	// 4
	mod_idx = 4;
	intt_keyswitch_1_para1(target[0], target_INTT_one_p1, idt[0], sidt[0], modulus[0]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[1], target_INTT_one_p0, idt[1], sidt[1], modulus[1]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[2], target_INTT_one_p1, idt[2], sidt[2], modulus[2]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[3], target_INTT_one_p0, idt[3], sidt[3], modulus[3]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[4], target_INTT_one_p1, idt[4], sidt[4], modulus[4]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);
}

*/

/*
// 其实这个函数也不�??要用这么多资�??
// 片外�??要存原本密文，和临时存的sum，每�??
template <unsigned modcount, unsigned ntt1_para, unsigned intt2_para, unsigned ntt2_para>
void switchkey_k4_layer3(UDTYPE encrypted_in[modcount][2][BRAMNUM][BRAMSIZE],
						 UDTYPE encrypted[2][BRAMNUM][BRAMSIZE],
						 UDTYPE target[modcount][BRAMNUM][BRAMSIZE],
						 UDTYPEin *data_out0, UDTYPEin *keys_in0, UDTYPEin *sum_in0, UDTYPEin *out_0,
						 UDTYPE rp[ROOTMODCOUNT][N],
						 UDTYPE srp[ROOTMODCOUNT][N],
						 UDTYPE idt[ROOTMODCOUNT][N],
						 UDTYPE sidt[ROOTMODCOUNT][N],
						 UDTYPE modulus[ROOTMODCOUNT],
						 UDTYPE modulus_ratio[ROOTMODCOUNT][2],
						 UDTYPE modswitch_factors[modcount],
						 uint32_t galois_elt[12],
						 UDTYPE target_INTT_one_p0[BRAMNUM][BRAMSIZE],
						 UDTYPE target_INTT_one_p1[BRAMNUM][BRAMSIZE],
						 UDTYPE buffer_NTT0_one_p0[BRAMNUM][BRAMSIZE],
						 UDTYPE buffer_NTT0_one_p1[BRAMNUM][BRAMSIZE],
						 UDTYPE buffer_NTT1_p0[2][BRAMNUM][BRAMSIZE],
						 UDTYPE buffer_NTT1_p1[2][BRAMNUM][BRAMSIZE],
						 UDTYPE buffer_dyadic_last[2][BRAMNUM][BRAMSIZE],
						 UDTYPE buffer_dyadic_p0[2][BRAMNUM][BRAMSIZE],
						 UDTYPE buffer_dyadic_p1[2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE
	UDTYPE coeff_count = N;
	UDTYPE decomp_mod_count = modcount;
	UDTYPE rns_mod_count = modcount + 1;
	UDTYPE key_mod_count = ROOTMODCOUNT;
	UDTYPE modulus_all[ROOTMODCOUNT];

	UDTYPE half = modulus[INITMODCOUNT - 1] >> 1;

	UDTYPE key_index = 1;

	for (int cout = 0; cout < 3; cout++)
	{
#pragma HLS UNROLL
		int mod_idx = 4;

		// m0
		galois_mc1(encrypted_in[0], encrypted, target[0], galois_elt[0]);
		data_out_mc1(encrypted, data_out0);

		intt_keyswitch_1_para1(target[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);

		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);

		key_index += 2 * N;

		// m1
		galois_mc1(encrypted_in[1], encrypted, target[1], galois_elt[0]);
		data_out_mc1(encrypted, data_out0);

		intt_keyswitch_1_para1(target[1], target_INTT_one_p1, idt[1], sidt[1], modulus[1]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;

		// m2
		galois_mc1(encrypted_in[2], encrypted, target[2], galois_elt[0]);
		data_out_mc1(encrypted, data_out0);

		intt_keyswitch_1_para1(target[2], target_INTT_one_p0, idt[2], sidt[2], modulus[2]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[5][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;

		// m3
		galois_mc1(encrypted_in[3], encrypted, target[3], galois_elt[0]);
		data_out_mc1(encrypted, data_out0);

		intt_keyswitch_1_para1(target[3], target_INTT_one_p1, idt[3], sidt[3], modulus[3]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;

		// intt
		intt_keyswitch_2<intt2_para>(
			buffer_dyadic_last, idt[mod_idx], sidt[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		// real m0
		mod_idx = 0;

		// m0
		intt_keyswitch_1_para1(target[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][mod_idx], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m1
		intt_keyswitch_1_para1(target[1], target_INTT_one_p1, idt[1], sidt[1], modulus[1]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m2
		intt_keyswitch_1_para1(target[2], target_INTT_one_p0, idt[2], sidt[2], modulus[2]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m3
		intt_keyswitch_1_para1(target[3], target_INTT_one_p1, idt[3], sidt[3], modulus[3]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(
			buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

		multiply_sub_add_cc_mc1(data_out0, sum_in0, out_0,
								buffer_NTT1_p0, buffer_dyadic_p0,
								modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

		// 1

		mod_idx = 1;
		// m0
		intt_keyswitch_1_para1(target[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m1
		intt_keyswitch_1_para1(target[1], target_INTT_one_p1, idt[1], sidt[1], modulus[1]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m2
		intt_keyswitch_1_para1(target[2], target_INTT_one_p0, idt[2], sidt[2], modulus[2]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m3
		intt_keyswitch_1_para1(target[3], target_INTT_one_p1, idt[3], sidt[3], modulus[3]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(
			buffer_dyadic_last, buffer_NTT1_p1, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

		multiply_sub_add_cc_mc1(data_out0, sum_in0, out_0,
								buffer_NTT1_p0, buffer_dyadic_p1,
								modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);
		//		multiply_sub<1, BRAMNUM, BRAMSIZE>(buffer_NTT1_p1, buffer_dyadic_p1,
		//				encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

		// 2
		mod_idx = 2;
		// m0
		intt_keyswitch_1_para1(target[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m1
		intt_keyswitch_1_para1(target[1], target_INTT_one_p1, idt[1], sidt[1], modulus[1]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m2
		intt_keyswitch_1_para1(target[2], target_INTT_one_p0, idt[2], sidt[2], modulus[2]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m3
		intt_keyswitch_1_para1(target[3], target_INTT_one_p1, idt[3], sidt[3], modulus[3]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(
			buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

		multiply_sub_add_cc_mc1(data_out0, sum_in0, out_0,
								buffer_NTT1_p0, buffer_dyadic_p0,
								modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);
		//		multiply_sub<1, BRAMNUM, BRAMSIZE>(buffer_NTT1_p0, buffer_dyadic_p0,
		//				encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

		// 3

		mod_idx = 3;
		// m0
		intt_keyswitch_1_para1(target[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m1
		intt_keyswitch_1_para1(target[1], target_INTT_one_p1, idt[1], sidt[1], modulus[1]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m2
		intt_keyswitch_1_para1(target[2], target_INTT_one_p0, idt[2], sidt[2], modulus[2]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;
		// m3
		intt_keyswitch_1_para1(target[3], target_INTT_one_p1, idt[3], sidt[3], modulus[3]);
		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(
			buffer_dyadic_last, buffer_NTT1_p1, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

		multiply_sub_add_cc_mc1(data_out0, sum_in0, out_0,
								buffer_NTT1_p0, buffer_dyadic_p1,
								modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);
		//		multiply_sub<1, BRAMNUM, BRAMSIZE>(buffer_NTT1_p1, buffer_dyadic_p1,
		//				 encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);
	}
}
*/

template <unsigned modcount, unsigned ntt1_para, unsigned intt2_para, unsigned ntt2_para>
void switchkey_k3_new(
	// UDTYPE encrypted_in[modcount][2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted0[modcount][2][BRAMNUM][BRAMSIZE],
	UDTYPE target0[modcount][BRAMNUM][BRAMSIZE],
	// UDTYPEin *data_in,
	UDTYPEin *keys_in0, UDTYPEin *data_out,
	UDTYPE rp[INITMODCOUNT][RPBRAMNUM][RPBRAMSIZE],
	UDTYPE srp[INITMODCOUNT][RPBRAMNUM][RPBRAMSIZE],
	UDTYPE idt[INITMODCOUNT][RPBRAMNUM][RPBRAMSIZE],
	UDTYPE sidt[INITMODCOUNT][RPBRAMNUM][RPBRAMSIZE],
	UDTYPE modulus[INITMODCOUNT],
	UDTYPE modcopy[INITMODCOUNT],
	UDTYPE modulus_ratio[INITMODCOUNT][3],
	UDTYPE modswitch_factors[modcount],
	uint32_t galois_elt[12],
	UDTYPE target_INTT[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT0[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT0_[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT1_p0[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT1_p1[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_last[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_p0[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_p1[2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	UDTYPE half = modulus[INITMODCOUNT - 1] >> 1;

	UDTYPE key_index = 0;
	UDTYPE data_index = 0;
	size_t index, j;

	// c <= 24
//	 data_in_mc3(data_in + data_index, encrypted_in);
	// galois<modcount>(encrypted_in, encrypted0, target0, galois_elt[0]);
	// data_index += 2 * N;

	// j = modcount
	j = modcount;									// j = 3
	index = (j == modcount ? INITMODCOUNT - 1 : j); // index = 3

	intt_keyswitch_1_new<modcount>(target0, target_INTT, idt, sidt, modulus);

	ntt_keyswitch_1_new<modcount, ntt1_para>(
		target_INTT, buffer_NTT0, rp[index], srp[index], modcopy, modulus_ratio[index][0], modulus_ratio[index][2]);

	// target0[2]只为了占�?
	dyadic_unpara_poly_new<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_last, modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], j, target0[2]);

	key_index += 2 * N;

	intt_keyswitch_2_new<intt2_para>(
		buffer_dyadic_last, idt[index], sidt[index], modulus_ratio[index][0], modulus_ratio[index][2]);

	// 0
	j = 0;											// j = 0
	index = (j == modcount ? INITMODCOUNT - 1 : j); // index = 0
	ntt_keyswitch_1_new<modcount, ntt1_para>(
		target_INTT, buffer_NTT0_, rp[index], srp[index], modcopy, modulus_ratio[index][0], modulus_ratio[index][2]);

	dyadic_unpara_poly_new<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p0, modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], j, target0[j]);

	key_index += 2 * N;

	ntt_keyswitch_2_new<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
								   rp[index], srp[index], modulus_ratio[index][0], modulus_ratio[index][2], half);


	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted0[index], modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], modswitch_factors[index]);


	// 1
	j = 1;											// j = 1
	index = (j == modcount ? INITMODCOUNT - 1 : j); // index = 1
	ntt_keyswitch_1_new<modcount, ntt1_para>(
		target_INTT, buffer_NTT0, rp[index], srp[index], modcopy, modulus_ratio[index][0], modulus_ratio[index][2]);

	dyadic_unpara_poly_new<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_p1, modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], j, target0[j]);

	key_index += 2 * N;

	ntt_keyswitch_2_new<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
								   rp[index], srp[index], modulus_ratio[index][0], modulus_ratio[index][2], half);

	multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
					encrypted0[index], modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], modswitch_factors[index]);

	// 2
	j = 2;											// j = 2
	index = (j == modcount ? INITMODCOUNT - 1 : j); // index = 2
	ntt_keyswitch_1_new<modcount, ntt1_para>(
		target_INTT, buffer_NTT0_, rp[index], srp[index], modcopy, modulus_ratio[index][0], modulus_ratio[index][2]);

	dyadic_unpara_poly_new<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p0, modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], j, target0[j]);

	key_index += 2 * N;

	ntt_keyswitch_2_new<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
								   rp[index], srp[index], modulus_ratio[index][0], modulus_ratio[index][2], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted0[index], modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], modswitch_factors[index]);


	j = 3;											// j = 2
	index = (j == modcount ? INITMODCOUNT - 1 : j); // index = 2
	ntt_keyswitch_1_new<modcount, ntt1_para>(
		target_INTT, buffer_NTT0, rp[index], srp[index], modcopy, modulus_ratio[index][0], modulus_ratio[index][2]);

	dyadic_unpara_poly_new<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_p1, modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], j, target0[j]);

	key_index += 2 * N;

	ntt_keyswitch_2_new<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
								   rp[index], srp[index], modulus_ratio[index][0], modulus_ratio[index][2], half);

	multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
					encrypted0[index], modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], modswitch_factors[index]);
	data_out_mc3(encrypted0, data_out);
	data_out += 2 * N;
}

template <unsigned modcount, unsigned ntt1_para, unsigned intt2_para, unsigned ntt2_para>
void rotate_k3_new(
	UDTYPE encrypted_in[modcount][2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted0[modcount][2][BRAMNUM][BRAMSIZE],
	UDTYPE target0[modcount][BRAMNUM][BRAMSIZE],
	UDTYPEin *data_in,
	UDTYPEin *keys_in0, UDTYPEin *data_out,
	UDTYPE rp[INITMODCOUNT][RPBRAMNUM][RPBRAMSIZE],
	UDTYPE srp[INITMODCOUNT][RPBRAMNUM][RPBRAMSIZE],
	UDTYPE idt[INITMODCOUNT][RPBRAMNUM][RPBRAMSIZE],
	UDTYPE sidt[INITMODCOUNT][RPBRAMNUM][RPBRAMSIZE],
	UDTYPE modulus[INITMODCOUNT],
	UDTYPE modcopy[INITMODCOUNT],
	UDTYPE modulus_ratio[INITMODCOUNT][3],
	UDTYPE modswitch_factors[modcount],
	uint32_t galois_elt[12],
	UDTYPE target_INTT[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT0[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT0_[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT1_p0[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT1_p1[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_last[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_p0[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_p1[2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	UDTYPE half = modulus[INITMODCOUNT - 1] >> 1;

	UDTYPE key_index = 0;
	UDTYPE data_index = 0;
	size_t index, j;

	// c <= 24
	 //data_in_mc4(data_in + data_index, encrypted_in);
	 galois<modcount>(encrypted_in, encrypted0, target0, galois_elt[0]);
	 data_index += 2 * N;

	// j = modcount
	j = modcount;									// j = 3
	index = (j == modcount ? INITMODCOUNT - 1 : j); // index = 3

	intt_keyswitch_1_new<modcount>(target0, target_INTT, idt, sidt, modulus);

	ntt_keyswitch_1_new<modcount, ntt1_para>(
		target_INTT, buffer_NTT0, rp[index], srp[index], modcopy, modulus_ratio[index][0], modulus_ratio[index][2]);

	// target0[2]只为了占�?
	dyadic_unpara_poly_new<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_last, modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], j, target0[2]);

	key_index += 2 * N;

	intt_keyswitch_2_new<intt2_para>(
		buffer_dyadic_last, idt[index], sidt[index], modulus_ratio[index][0], modulus_ratio[index][2]);

	// 0
	j = 0;											// j = 0
	index = (j == modcount ? INITMODCOUNT - 1 : j); // index = 0
	ntt_keyswitch_1_new<modcount, ntt1_para>(
		target_INTT, buffer_NTT0_, rp[index], srp[index], modcopy, modulus_ratio[index][0], modulus_ratio[index][2]);

	dyadic_unpara_poly_new<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p0, modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], j, target0[j]);

	key_index += 2 * N;

	ntt_keyswitch_2_new<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
								   rp[index], srp[index], modulus_ratio[index][0], modulus_ratio[index][2], half);


	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted0[index], modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], modswitch_factors[index]);


	// 1
	j = 1;											// j = 1
	index = (j == modcount ? INITMODCOUNT - 1 : j); // index = 1
	ntt_keyswitch_1_new<modcount, ntt1_para>(
		target_INTT, buffer_NTT0, rp[index], srp[index], modcopy, modulus_ratio[index][0], modulus_ratio[index][2]);

	dyadic_unpara_poly_new<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_p1, modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], j, target0[j]);

	key_index += 2 * N;

	ntt_keyswitch_2_new<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
								   rp[index], srp[index], modulus_ratio[index][0], modulus_ratio[index][2], half);

	multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
					encrypted0[index], modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], modswitch_factors[index]);

	// 2
	j = 2;											// j = 2
	index = (j == modcount ? INITMODCOUNT - 1 : j); // index = 2
	ntt_keyswitch_1_new<modcount, ntt1_para>(
		target_INTT, buffer_NTT0_, rp[index], srp[index], modcopy, modulus_ratio[index][0], modulus_ratio[index][2]);

	dyadic_unpara_poly_new<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p0, modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], j, target0[j]);

	key_index += 2 * N;

	ntt_keyswitch_2_new<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
								   rp[index], srp[index], modulus_ratio[index][0], modulus_ratio[index][2], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted0[index], modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], modswitch_factors[index]);


	j = 3;											// j = 2
	index = (j == modcount ? INITMODCOUNT - 1 : j); // index = 2
	ntt_keyswitch_1_new<modcount, ntt1_para>(
		target_INTT, buffer_NTT0, rp[index], srp[index], modcopy, modulus_ratio[index][0], modulus_ratio[index][2]);

	dyadic_unpara_poly_new<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_p1, modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], j, target0[j]);

	key_index += 2 * N;

	ntt_keyswitch_2_new<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
								   rp[index], srp[index], modulus_ratio[index][0], modulus_ratio[index][2], half);

	multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
					encrypted0[index], modulus_ratio[index][0], modulus_ratio[index][1], modulus_ratio[index][2], modswitch_factors[index]);
//	data_out_mc3(encrypted0, data_out);
//	data_out += 2 * N;
}

/*
template <unsigned modcount, unsigned ntt1_para, unsigned intt2_para, unsigned ntt2_para>
void switchkey_k3_layer3(
	UDTYPE encrypted_in[modcount][2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted0[modcount][2][BRAMNUM][BRAMSIZE],
	UDTYPE target0[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted1[modcount][2][BRAMNUM][BRAMSIZE],
	UDTYPE target1[modcount][BRAMNUM][BRAMSIZE],
	UDTYPEin *data_in, UDTYPEin *keys_in0, UDTYPEin *data_out,
	UDTYPE rp[ROOTMODCOUNT][N],
	UDTYPE srp[ROOTMODCOUNT][N],
	UDTYPE idt[ROOTMODCOUNT][N],
	UDTYPE sidt[ROOTMODCOUNT][N],
	UDTYPE modulus[ROOTMODCOUNT],
	UDTYPE modulus_ratio[ROOTMODCOUNT][2],
	UDTYPE modswitch_factors[modcount],
	uint32_t galois_elt[12],
	UDTYPE target_INTT[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT0[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT0_[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT1_p0[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT1_p1[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_last[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_p0[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_p1[2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	UDTYPE half = modulus[INITMODCOUNT - 1] >> 1;

	UDTYPE key_index = 0;
	UDTYPE data_index = 0;

	// c <= 24
	for (int c = 0; c < 24; c += 2)
	{
#pragma HLS UNROLL

		data_in_mc3(data_in + data_index, encrypted_in);
		galois<modcount>(encrypted_in, encrypted0, target0, galois_elt[0]);
		data_index += 2 * N;

		// before
		intt_keyswitch_1<modcount>(target0, target_INTT, idt, sidt, modulus);

		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0, rp[3], srp[3], modulus[3], modulus_ratio[3][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_last, modulus[3], modulus_ratio[3][0], modulus_ratio[3][1]);

		key_index += 2 * N;

		intt_keyswitch_2<intt2_para>(
			buffer_dyadic_last, idt[3], sidt[3], modulus[3], modulus_ratio[3][1]);

		// 0
		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0_, rp[0], srp[0], modulus[0], modulus_ratio[0][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p0, modulus[0], modulus_ratio[0][0], modulus_ratio[0][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
								   rp[0], srp[0], modulus[0], modulus_ratio[0][1], half);

		multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
						encrypted0[0], modulus[0], modulus_ratio[0][0], modulus_ratio[0][1], modswitch_factors[0]);

		// 1
		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0, rp[1], srp[1], modulus[1], modulus_ratio[1][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_p1, modulus[1], modulus_ratio[1][0], modulus_ratio[1][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
								   rp[1], srp[1], modulus[1], modulus_ratio[1][1], half);

		multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
						encrypted0[1], modulus[1], modulus_ratio[1][0], modulus_ratio[1][1], modswitch_factors[1]);

		// 2

		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0_, rp[2], srp[2], modulus[2], modulus_ratio[2][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p0, modulus[2], modulus_ratio[2][0], modulus_ratio[2][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
								   rp[2], srp[2], modulus[2], modulus_ratio[2][1], half);

		multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
						encrypted0[2], modulus[2], modulus_ratio[2][0], modulus_ratio[2][1], modswitch_factors[2]);

		data_out_mc3(encrypted0, data_out);
		data_out += 2 * N;

		// rotation 2 -----------------------------------------------------------------------

		data_in_mc3(data_in + data_index, encrypted_in);
		galois<modcount>(encrypted_in, encrypted1, target1, galois_elt[0]);
		data_index += 2 * N;

		// before
		intt_keyswitch_1<modcount>(target1, target_INTT, idt, sidt, modulus);

		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0, rp[3], srp[3], modulus[3], modulus_ratio[3][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_last, modulus[3], modulus_ratio[3][0], modulus_ratio[3][1]);

		key_index += 2 * N;

		intt_keyswitch_2<intt2_para>(
			buffer_dyadic_last, idt[3], sidt[3], modulus[3], modulus_ratio[3][1]);

		// 0
		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0_, rp[0], srp[0], modulus[0], modulus_ratio[0][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p1, modulus[0], modulus_ratio[0][0], modulus_ratio[0][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
								   rp[0], srp[0], modulus[0], modulus_ratio[0][1], half);

		multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
						encrypted1[0], modulus[0], modulus_ratio[0][0], modulus_ratio[0][1], modswitch_factors[0]);

		// 1
		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0, rp[1], srp[1], modulus[1], modulus_ratio[1][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_p0, modulus[1], modulus_ratio[1][0], modulus_ratio[1][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
								   rp[1], srp[1], modulus[1], modulus_ratio[1][1], half);

		multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
						encrypted1[1], modulus[1], modulus_ratio[1][0], modulus_ratio[1][1], modswitch_factors[1]);

		// 2

		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0_, rp[2], srp[2], modulus[2], modulus_ratio[2][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p1, modulus[2], modulus_ratio[2][0], modulus_ratio[2][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
								   rp[2], srp[2], modulus[2], modulus_ratio[2][1], half);

		multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
						encrypted1[2], modulus[2], modulus_ratio[2][0], modulus_ratio[2][1], modswitch_factors[2]);

		data_out_mc3(encrypted1, data_out);
		data_out += 2 * N;
	}

	// the 25 rotation
	data_in_mc3(data_in + data_index, encrypted_in);
	galois<modcount>(encrypted_in, encrypted0, target0, galois_elt[0]);
	data_index += 2 * N;

	// before
	intt_keyswitch_1<modcount>(target0, target_INTT, idt, sidt, modulus);

	ntt_keyswitch_1<modcount, ntt1_para>(
		target_INTT, buffer_NTT0, rp[3], srp[3], modulus[3], modulus_ratio[3][1]);

	dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
								 buffer_dyadic_last, modulus[3], modulus_ratio[3][0], modulus_ratio[3][1]);

	key_index += 2 * N;

	intt_keyswitch_2<intt2_para>(
		buffer_dyadic_last, idt[3], sidt[3], modulus[3], modulus_ratio[3][1]);

	// 0
	ntt_keyswitch_1<modcount, ntt1_para>(
		target_INTT, buffer_NTT0_, rp[0], srp[0], modulus[0], modulus_ratio[0][1]);

	dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
								 buffer_dyadic_p0, modulus[0], modulus_ratio[0][0], modulus_ratio[0][1]);

	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
							   rp[0], srp[0], modulus[0], modulus_ratio[0][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted0[0], modulus[0], modulus_ratio[0][0], modulus_ratio[0][1], modswitch_factors[0]);

	// 1
	ntt_keyswitch_1<modcount, ntt1_para>(
		target_INTT, buffer_NTT0, rp[1], srp[1], modulus[1], modulus_ratio[1][1]);

	dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
								 buffer_dyadic_p1, modulus[1], modulus_ratio[1][0], modulus_ratio[1][1]);

	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
							   rp[1], srp[1], modulus[1], modulus_ratio[1][1], half);

	multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
					encrypted0[1], modulus[1], modulus_ratio[1][0], modulus_ratio[1][1], modswitch_factors[1]);

	// 2

	ntt_keyswitch_1<modcount, ntt1_para>(
		target_INTT, buffer_NTT0_, rp[2], srp[2], modulus[2], modulus_ratio[2][1]);

	dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
								 buffer_dyadic_p0, modulus[2], modulus_ratio[2][0], modulus_ratio[2][1]);

	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
							   rp[2], srp[2], modulus[2], modulus_ratio[2][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted0[2], modulus[2], modulus_ratio[2][0], modulus_ratio[2][1], modswitch_factors[2]);

	data_out_mc3(encrypted0, data_out);
	data_out += 2 * N;
}


template <unsigned modcount, unsigned ntt1_para, unsigned intt2_para, unsigned ntt2_para>
void switchkey_k3_layer3_copy(
	UDTYPE encrypted_in[modcount][2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted0[modcount][2][BRAMNUM][BRAMSIZE],
	UDTYPE target0[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted1[modcount][2][BRAMNUM][BRAMSIZE],
	UDTYPE target1[modcount][BRAMNUM][BRAMSIZE],
	UDTYPEin *data_in, UDTYPEin *keys_in0, UDTYPEin *data_out,
	UDTYPE rp[ROOTMODCOUNT][N],
	UDTYPE srp[ROOTMODCOUNT][N],
	UDTYPE idt[ROOTMODCOUNT][N],
	UDTYPE sidt[ROOTMODCOUNT][N],
	UDTYPE modulus[ROOTMODCOUNT],
	UDTYPE modulus_ratio[ROOTMODCOUNT][2],
	UDTYPE modswitch_factors[modcount],
	uint32_t galois_elt[12],
	UDTYPE target_INTT[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT0[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT0_[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT1_p0[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT1_p1[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_last[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_p0[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_p1[2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	UDTYPE half = modulus[INITMODCOUNT - 1] >> 1;

	UDTYPE key_index = 0;
	UDTYPE data_index = 0;

	// c <= 24
	for (int c = 0; c < 24; c += 2)
	{
#pragma HLS UNROLL

		data_in_mc3(data_in + data_index, encrypted_in);
		galois<modcount>(encrypted_in, encrypted0, target0, galois_elt[0]);
		data_index += 2 * N;

		// before
		intt_keyswitch_1<modcount>(target0, target_INTT, idt, sidt, modulus);

		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0, rp[3], srp[3], modulus[3], modulus_ratio[3][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_last, modulus[3], modulus_ratio[3][0], modulus_ratio[3][1]);

		key_index += 2 * N;

		intt_keyswitch_2<intt2_para>(
			buffer_dyadic_last, idt[3], sidt[3], modulus[3], modulus_ratio[3][1]);

		// 0
		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0_, rp[0], srp[0], modulus[0], modulus_ratio[0][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p0, modulus[0], modulus_ratio[0][0], modulus_ratio[0][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
								   rp[0], srp[0], modulus[0], modulus_ratio[0][1], half);

		multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
						encrypted0[0], modulus[0], modulus_ratio[0][0], modulus_ratio[0][1], modswitch_factors[0]);

		// 1
		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0, rp[1], srp[1], modulus[1], modulus_ratio[1][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_p1, modulus[1], modulus_ratio[1][0], modulus_ratio[1][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
								   rp[1], srp[1], modulus[1], modulus_ratio[1][1], half);

		multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
						encrypted0[1], modulus[1], modulus_ratio[1][0], modulus_ratio[1][1], modswitch_factors[1]);

		// 2

		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0_, rp[2], srp[2], modulus[2], modulus_ratio[2][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p0, modulus[2], modulus_ratio[2][0], modulus_ratio[2][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
								   rp[2], srp[2], modulus[2], modulus_ratio[2][1], half);

		multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
						encrypted0[2], modulus[2], modulus_ratio[2][0], modulus_ratio[2][1], modswitch_factors[2]);

		data_out_mc3(encrypted0, data_out);
		data_out += 2 * N;

		// rotation 2 -----------------------------------------------------------------------

		data_in_mc3(data_in + data_index, encrypted_in);
		galois<modcount>(encrypted_in, encrypted1, target1, galois_elt[0]);
		data_index += 2 * N;

		// before
		intt_keyswitch_1<modcount>(target1, target_INTT, idt, sidt, modulus);

		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0, rp[3], srp[3], modulus[3], modulus_ratio[3][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_last, modulus[3], modulus_ratio[3][0], modulus_ratio[3][1]);

		key_index += 2 * N;

		intt_keyswitch_2<intt2_para>(
			buffer_dyadic_last, idt[3], sidt[3], modulus[3], modulus_ratio[3][1]);

		// 0
		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0_, rp[0], srp[0], modulus[0], modulus_ratio[0][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p1, modulus[0], modulus_ratio[0][0], modulus_ratio[0][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
								   rp[0], srp[0], modulus[0], modulus_ratio[0][1], half);

		multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
						encrypted1[0], modulus[0], modulus_ratio[0][0], modulus_ratio[0][1], modswitch_factors[0]);

		// 1
		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0, rp[1], srp[1], modulus[1], modulus_ratio[1][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
									 buffer_dyadic_p0, modulus[1], modulus_ratio[1][0], modulus_ratio[1][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
								   rp[1], srp[1], modulus[1], modulus_ratio[1][1], half);

		multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
						encrypted1[1], modulus[1], modulus_ratio[1][0], modulus_ratio[1][1], modswitch_factors[1]);

		// 2

		ntt_keyswitch_1<modcount, ntt1_para>(
			target_INTT, buffer_NTT0_, rp[2], srp[2], modulus[2], modulus_ratio[2][1]);

		dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
									 buffer_dyadic_p1, modulus[2], modulus_ratio[2][0], modulus_ratio[2][1]);

		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
								   rp[2], srp[2], modulus[2], modulus_ratio[2][1], half);

		multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
						encrypted1[2], modulus[2], modulus_ratio[2][0], modulus_ratio[2][1], modswitch_factors[2]);

		data_out_mc3(encrypted1, data_out);
		data_out += 2 * N;
	}

	// the 25 rotation
	data_in_mc3(data_in + data_index, encrypted_in);
	galois<modcount>(encrypted_in, encrypted0, target0, galois_elt[0]);
	data_index += 2 * N;

	// before
	intt_keyswitch_1<modcount>(target0, target_INTT, idt, sidt, modulus);

	ntt_keyswitch_1<modcount, ntt1_para>(
		target_INTT, buffer_NTT0, rp[3], srp[3], modulus[3], modulus_ratio[3][1]);

	dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
								 buffer_dyadic_last, modulus[3], modulus_ratio[3][0], modulus_ratio[3][1]);

	key_index += 2 * N;

	intt_keyswitch_2<intt2_para>(
		buffer_dyadic_last, idt[3], sidt[3], modulus[3], modulus_ratio[3][1]);

	// 0
	ntt_keyswitch_1<modcount, ntt1_para>(
		target_INTT, buffer_NTT0_, rp[0], srp[0], modulus[0], modulus_ratio[0][1]);

	dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
								 buffer_dyadic_p0, modulus[0], modulus_ratio[0][0], modulus_ratio[0][1]);

	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
							   rp[0], srp[0], modulus[0], modulus_ratio[0][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted0[0], modulus[0], modulus_ratio[0][0], modulus_ratio[0][1], modswitch_factors[0]);

	// 1
	ntt_keyswitch_1<modcount, ntt1_para>(
		target_INTT, buffer_NTT0, rp[1], srp[1], modulus[1], modulus_ratio[1][1]);

	dyadic_unpara_poly<modcount>(buffer_NTT0, keys_in0 + key_index,
								 buffer_dyadic_p1, modulus[1], modulus_ratio[1][0], modulus_ratio[1][1]);

	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p1,
							   rp[1], srp[1], modulus[1], modulus_ratio[1][1], half);

	multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
					encrypted0[1], modulus[1], modulus_ratio[1][0], modulus_ratio[1][1], modswitch_factors[1]);

	// 2

	ntt_keyswitch_1<modcount, ntt1_para>(
		target_INTT, buffer_NTT0_, rp[2], srp[2], modulus[2], modulus_ratio[2][1]);

	dyadic_unpara_poly<modcount>(buffer_NTT0_, keys_in0 + key_index,
								 buffer_dyadic_p0, modulus[2], modulus_ratio[2][0], modulus_ratio[2][1]);

	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(buffer_dyadic_last, buffer_NTT1_p0,
							   rp[2], srp[2], modulus[2], modulus_ratio[2][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted0[2], modulus[2], modulus_ratio[2][0], modulus_ratio[2][1], modswitch_factors[2]);

	data_out_mc3(encrypted0, data_out);
	data_out += 2 * N;
}

*/

/*
template <unsigned modcount, unsigned ntt1_para, unsigned intt2_para, unsigned ntt2_para>
void switchkey_k3_para1(UDTYPE encrypted[modcount][2][BRAMNUM][BRAMSIZE],
						UDTYPE target[modcount][BRAMNUM][BRAMSIZE],
						UDTYPEin *keys_in0,
						UDTYPE rp[ROOTMODCOUNT][N],
						UDTYPE srp[ROOTMODCOUNT][N],
						UDTYPE idt[ROOTMODCOUNT][N],
						UDTYPE sidt[ROOTMODCOUNT][N],
						UDTYPE modulus[ROOTMODCOUNT],
						UDTYPE modulus_ratio[ROOTMODCOUNT][2],
						UDTYPE modswitch_factors[modcount],
						UDTYPE target_INTT_one_p0[BRAMNUM][BRAMSIZE],
						UDTYPE target_INTT_one_p1[BRAMNUM][BRAMSIZE],
						UDTYPE buffer_NTT0_one_p0[BRAMNUM][BRAMSIZE],
						UDTYPE buffer_NTT0_one_p1[BRAMNUM][BRAMSIZE],
						UDTYPE buffer_NTT1_p0[2][BRAMNUM][BRAMSIZE],
						UDTYPE buffer_NTT1_p1[2][BRAMNUM][BRAMSIZE],
						UDTYPE buffer_dyadic_last[2][BRAMNUM][BRAMSIZE],
						UDTYPE buffer_dyadic_p0[2][BRAMNUM][BRAMSIZE],
						UDTYPE buffer_dyadic_p1[2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE

	UDTYPE half = modulus[INITMODCOUNT - 1] >> 1;

	UDTYPE key_index = 0;

	int mod_idx = 3;

	intt_keyswitch_1_para1(target[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);

	key_index += 2 * N;

	intt_keyswitch_1_para1(target[1], target_INTT_one_p1, idt[1], sidt[1], modulus[1]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[2], target_INTT_one_p0, idt[2], sidt[2], modulus[2]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[5][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_2<intt2_para>(
		buffer_dyadic_last, idt[mod_idx], sidt[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	// 0

	mod_idx = 0;

	intt_keyswitch_1_para1(target[0], target_INTT_one_p1, idt[0], sidt[0], modulus[0]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][mod_idx], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[1], target_INTT_one_p0, idt[1], sidt[1], modulus[1]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[2], target_INTT_one_p1, idt[2], sidt[2], modulus[2]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

	// 1

	mod_idx = 1;

	intt_keyswitch_1_para1(target[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[1], target_INTT_one_p1, idt[1], sidt[1], modulus[1]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[2], target_INTT_one_p0, idt[2], sidt[2], modulus[2]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p1, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p1, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p1, buffer_dyadic_p1,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);

	// 2
	mod_idx = 2;
	intt_keyswitch_1_para1(target[0], target_INTT_one_p1, idt[0], sidt[0], modulus[0]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[1], target_INTT_one_p0, idt[1], sidt[1], modulus[1]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	intt_keyswitch_1_para1(target[2], target_INTT_one_p1, idt[2], sidt[2], modulus[2]);
	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p1, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1], modswitch_factors[mod_idx]);
}

*/

/*
template <unsigned modcount, unsigned ntt1_para, unsigned intt2_para, unsigned ntt2_para>
void switchkey_k1_layer5(
	UDTYPE encrypted_in[2][BRAMNUM][BRAMSIZE], //��
	UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],   //��
	UDTYPE target0[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE], //��
	UDTYPE target1[modcount][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted2[2][BRAMNUM][BRAMSIZE],
	UDTYPEin *data_in, UDTYPEin *keys_in0, UDTYPEin *data_out,
	UDTYPE rp[ROOTMODCOUNT][N],
	UDTYPE srp[ROOTMODCOUNT][N],
	UDTYPE idt[ROOTMODCOUNT][N],
	UDTYPE sidt[ROOTMODCOUNT][N],
	UDTYPE modulus[ROOTMODCOUNT],
	UDTYPE modulus_ratio[ROOTMODCOUNT][2],
	UDTYPE modswitch_factors[modcount],
	uint32_t galois_elt[12],
	UDTYPE target_INTT_one_p0[BRAMNUM][BRAMSIZE],
	// UDTYPE target_INTT_one_p1[BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT0_one_p0[BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT0_one_p1[BRAMNUM][BRAMSIZE],
	UDTYPE buffer_NTT1_p0[2][BRAMNUM][BRAMSIZE],
	// UDTYPE buffer_NTT1_p1[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_last[2][BRAMNUM][BRAMSIZE],
	UDTYPE buffer_dyadic_p0[2][BRAMNUM][BRAMSIZE]
	// UDTYPE buffer_dyadic_p1[2][BRAMNUM][BRAMSIZE]
)
{
#pragma HLS INLINE

	UDTYPE half = modulus[INITMODCOUNT - 1] >> 1;

	UDTYPE key_index = 0;
	UDTYPE data_index = 0;

	// c <= 10
	for (int c = 0; c < 9; c += 3)
	{
#pragma HLS UNROLL

		int mod_idx = 1;

		data_in_mc1(data_in + data_index, encrypted_in);
		galois_mc1(encrypted_in, encrypted0, target0[0], galois_elt[0]); //������
		data_index += 2 * N;

		// before
		intt_keyswitch_1_para1(target0[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);

		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);

		key_index += 2 * N;

		intt_keyswitch_2<intt2_para>(
			buffer_dyadic_last, idt[mod_idx], sidt[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		// 0
		mod_idx = 0;

		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][mod_idx], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(
			buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

		multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
						encrypted0, modulus[0], modulus_ratio[0][0], modulus_ratio[0][1], modswitch_factors[0]); // encrypted0��

		data_out_mc1(encrypted0, data_out); //������
		data_out += 2 * N;

		// rotation 2_____________________________________________________________________

		mod_idx = 1;

		data_in_mc1(data_in + data_index, encrypted_in);
		galois_mc1(encrypted_in, encrypted1, target1[0], galois_elt[0]); //������
		data_index += 2 * N;

		// before
		intt_keyswitch_1_para1(target1[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);

		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);

		key_index += 2 * N;

		intt_keyswitch_2<intt2_para>(
			buffer_dyadic_last, idt[mod_idx], sidt[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		// 0
		mod_idx = 0;

		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][mod_idx], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(
			buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

		multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
						encrypted1, modulus[0], modulus_ratio[0][0], modulus_ratio[0][1], modswitch_factors[0]); // encrypted0��

		data_out_mc1(encrypted1, data_out); //������
		data_out += 2 * N;

		// rotation3 ________________________________________________________________________
		mod_idx = 1;

		data_in_mc1(data_in + data_index, encrypted_in);
		galois_mc1(encrypted_in, encrypted2, target0[0], galois_elt[0]); //������
		data_index += 2 * N;

		// before
		intt_keyswitch_1_para1(target0[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);

		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
				   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);

		key_index += 2 * N;

		intt_keyswitch_2<intt2_para>(
			buffer_dyadic_last, idt[mod_idx], sidt[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		// 0
		mod_idx = 0;

		ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p1,
										 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

		dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
				   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][mod_idx], modulus_ratio[mod_idx][1]);
		key_index += 2 * N;

		ntt_keyswitch_2<ntt2_para>(
			buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

		multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
						encrypted2, modulus[0], modulus_ratio[0][0], modulus_ratio[0][1], modswitch_factors[0]); // encrypted0��

		data_out_mc1(encrypted2, data_out); //������
		data_out += 2 * N;
	}
	// 10
	int mod_idx = 1;

	data_in_mc1(data_in + data_index, encrypted_in);
	galois_mc1(encrypted_in, encrypted0, target0[0], galois_elt[0]); //������
	data_index += 2 * N;

	// before
	intt_keyswitch_1_para1(target0[0], target_INTT_one_p0, idt[0], sidt[0], modulus[0]);

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p0,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p0, keys_in0 + key_index,
			   buffer_dyadic_last, modulus[mod_idx], modulus_ratio[mod_idx][0], modulus_ratio[mod_idx][1]);

	key_index += 2 * N;

	intt_keyswitch_2<intt2_para>(
		buffer_dyadic_last, idt[mod_idx], sidt[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	// 0
	mod_idx = 0;

	ntt_keyswitch_1_para1<ntt1_para>(target_INTT_one_p0, buffer_NTT0_one_p1,
									 rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1]);

	dyadic_mc1(buffer_NTT0_one_p1, keys_in0 + key_index,
			   buffer_dyadic_p0, modulus[mod_idx], modulus_ratio[mod_idx][mod_idx], modulus_ratio[mod_idx][1]);
	key_index += 2 * N;

	ntt_keyswitch_2<ntt2_para>(
		buffer_dyadic_last, buffer_NTT1_p0, rp[mod_idx], srp[mod_idx], modulus[mod_idx], modulus_ratio[mod_idx][1], half);

	multiply_sub<1>(buffer_NTT1_p0, buffer_dyadic_p0,
					encrypted0, modulus[0], modulus_ratio[0][0], modulus_ratio[0][1], modswitch_factors[0]); // encrypted0��

	data_out_mc1(encrypted0, data_out); //������
	data_out += 2 * N;
}

*/
