#pragma once
#include "define.h"
using namespace std;

void rescale_core(UDTYPE in, UDTYPE encrypted_last_mod, UDTYPE half_mod, UDTYPE &out,
				  UDTYPE modulus, UDTYPE modulus_ratio_0, UDTYPE modulus_ratio_1, UDTYPE inv_last);

void intt_mc3_new(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],
				  UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE], UDTYPE encrypted2[2][BRAMNUM][BRAMSIZE],
				  UDTYPE idt0[RPBRAMNUM][RPBRAMSIZE], UDTYPE idt1[RPBRAMNUM][RPBRAMSIZE], UDTYPE idt2[RPBRAMNUM][RPBRAMSIZE],
				  UDTYPE sidt0[RPBRAMNUM][RPBRAMSIZE], UDTYPE sidt1[RPBRAMNUM][RPBRAMSIZE], UDTYPE sidt2[RPBRAMNUM][RPBRAMSIZE],
				  UDTYPE modulus0, UDTYPE modulus1, UDTYPE modulus2);
void ntt_mc3_new(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],
				 UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE], UDTYPE encrypted2[2][BRAMNUM][BRAMSIZE],
				 UDTYPE rp0[RPBRAMNUM][RPBRAMSIZE], UDTYPE rp1[RPBRAMNUM][RPBRAMSIZE], UDTYPE rp2[RPBRAMNUM][RPBRAMSIZE],
				 UDTYPE srp0[RPBRAMNUM][RPBRAMSIZE], UDTYPE srp1[RPBRAMNUM][RPBRAMSIZE], UDTYPE srp2[RPBRAMNUM][RPBRAMSIZE],
				 UDTYPE modulus0, UDTYPE modulus1, UDTYPE modulus2);

void intt_mc3(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],
			  UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE], UDTYPE encrypted2[2][BRAMNUM][BRAMSIZE],
			  UDTYPE idt0[N], UDTYPE idt1[N], UDTYPE idt2[N],
			  UDTYPE sidt0[N], UDTYPE sidt1[N], UDTYPE sidt2[N],
			  UDTYPE modulus0, UDTYPE modulus1, UDTYPE modulus2);

void ntt_mc3(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],
			 UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE], UDTYPE encrypted2[2][BRAMNUM][BRAMSIZE],
			 UDTYPE rp0[N], UDTYPE rp1[N], UDTYPE rp2[N],
			 UDTYPE srp0[N], UDTYPE srp1[N], UDTYPE srp2[N],
			 UDTYPE modulus0, UDTYPE modulus1, UDTYPE modulus2);

void intt_mc2(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],
			  UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE],
			  UDTYPE idt0[N], UDTYPE idt1[N],
			  UDTYPE sidt0[N], UDTYPE sidt1[N],
			  UDTYPE modulus0, UDTYPE modulus1);

void ntt_mc2_new(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],
				 UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE],
				 UDTYPE rp0[RPBRAMNUM][RPBRAMSIZE], UDTYPE rp1[RPBRAMNUM][RPBRAMSIZE],
				 UDTYPE srp0[RPBRAMNUM][RPBRAMSIZE], UDTYPE srp1[RPBRAMNUM][RPBRAMSIZE],
				 UDTYPE modulus0, UDTYPE modulus1);

void ntt_mc2(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],
			 UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE],
			 UDTYPE rp0[N], UDTYPE rp1[N],
			 UDTYPE srp0[N], UDTYPE srp1[N],
			 UDTYPE modulus0, UDTYPE modulus1);

void ntt_mc1(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],
			 UDTYPE rp0[N],
			 UDTYPE srp0[N],
			 UDTYPE modulus0);

// BNC = BRAMNUM para count
template <unsigned BNC>
void rescale_para_5(UDTYPE encrypted_last[2][BRAMNUM][BRAMSIZE],
					UDTYPE encrypted_0[2][BRAMNUM][BRAMSIZE],
					UDTYPE encrypted_1[2][BRAMNUM][BRAMSIZE],
					UDTYPE encrypted_2[2][BRAMNUM][BRAMSIZE],
					UDTYPE encrypted_3[2][BRAMNUM][BRAMSIZE],
					UDTYPE encrypted_4[2][BRAMNUM][BRAMSIZE],
					UDTYPE half, UDTYPE half_mod[5],
					UDTYPE modulus_last, UDTYPE modulus_ratio_last[2],
					UDTYPE modulus[5], UDTYPE modulus_ratio[5][2], UDTYPE inv_last[5])
{

	UDTYPE last_temp;

	for (size_t j = 0; j < BRAMSIZE; j++)
	{
		for (size_t k = 0; k < 2; k++)
		{
			for (size_t i = 0; i < BRAMNUM; i += BNC)
			{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = encrypted_0 inter false
#pragma HLS DEPENDENCE variable = encrypted_1 inter false
#pragma HLS DEPENDENCE variable = encrypted_2 inter false
#pragma HLS DEPENDENCE variable = encrypted_3 inter false
#pragma HLS DEPENDENCE variable = encrypted_4 inter false
				for (size_t c = 0; c < BNC; c++)
				{

					size_t i_c;
					i_c = i + c;

					last_temp = barrett_reduce_63(encrypted_last[k][i_c][j] + half, modulus_last, modulus_ratio_last[1]);

					rescale_core(encrypted_0[k][i_c][j], last_temp, half_mod[0], encrypted_0[k][i_c][j],
								 modulus[0], modulus_ratio[0][0], modulus_ratio[0][1], inv_last[0]);
					rescale_core(encrypted_1[k][i_c][j], last_temp, half_mod[1], encrypted_1[k][i_c][j],
								 modulus[1], modulus_ratio[1][0], modulus_ratio[1][1], inv_last[1]);
					rescale_core(encrypted_2[k][i_c][j], last_temp, half_mod[2], encrypted_2[k][i_c][j],
								 modulus[2], modulus_ratio[2][0], modulus_ratio[2][1], inv_last[2]);
					rescale_core(encrypted_3[k][i_c][j], last_temp, half_mod[3], encrypted_3[k][i_c][j],
								 modulus[3], modulus_ratio[3][0], modulus_ratio[3][1], inv_last[3]);
					rescale_core(encrypted_4[k][i_c][j], last_temp, half_mod[4], encrypted_4[k][i_c][j],
								 modulus[4], modulus_ratio[4][0], modulus_ratio[4][1], inv_last[4]);
				}
			}
		}
	}
}

// BNC = BRAMNUM para count
template <unsigned BNC>
void rescale_para_mc3(UDTYPE encrypted_last[2][BRAMNUM][BRAMSIZE],
					  UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],
					  UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE],
					  UDTYPE encrypted2[2][BRAMNUM][BRAMSIZE],
					  UDTYPE half_last, UDTYPE half_mod0, UDTYPE half_mod1, UDTYPE half_mod2,
					  UDTYPE inv_last0, UDTYPE inv_last1, UDTYPE inv_last2,
					  UDTYPE modulus_fac_last[3],
					  UDTYPE modulus_fac0[3],
					  UDTYPE modulus_fac1[3],
					  UDTYPE modulus_fac2[3])
{

	UDTYPE last_temp;

	for (size_t j = 0; j < BRAMSIZE; j++)
	{
		for (size_t k = 0; k < 2; k++)
		{
			for (size_t i = 0; i < BRAMNUM; i += BNC)
			{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = encrypted0 inter false
#pragma HLS DEPENDENCE variable = encrypted1 inter false
#pragma HLS DEPENDENCE variable = encrypted2 inter false

				for (size_t c = 0; c < BNC; c++)
				{

					size_t i_c;
					i_c = i + c;

					// cout << i_c << endl;
					// cout << "encrypted_last[" << k << "][" << i_c << "][" << j << "] = " << encrypted_last[k][i_c][j] << endl;
					// cout << "half_last = " << half_last << endl;
					last_temp = barrett_reduce_63(encrypted_last[k][i_c][j] + half_last, modulus_fac_last[0], modulus_fac_last[2]);
					// cout << "last_temp = " << last_temp << endl;

					rescale_core(encrypted0[k][i_c][j], last_temp, half_mod0, encrypted0[k][i_c][j],
								 modulus_fac0[0], modulus_fac0[1], modulus_fac0[2], inv_last0);
					rescale_core(encrypted1[k][i_c][j], last_temp, half_mod1, encrypted1[k][i_c][j],
								 modulus_fac1[0], modulus_fac1[1], modulus_fac1[2], inv_last1);
					rescale_core(encrypted2[k][i_c][j], last_temp, half_mod2, encrypted2[k][i_c][j],
								 modulus_fac2[0], modulus_fac2[1], modulus_fac2[2], inv_last2);
				}
			}
		}
	}
}

// BNC = BRAMNUM para count
template <unsigned BNC>
void rescale_para_mc2(UDTYPE encrypted_last[2][BRAMNUM][BRAMSIZE],
					  UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],
					  UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE],
					  UDTYPE half_last, UDTYPE half_mod0, UDTYPE half_mod1,
					  UDTYPE inv_last0, UDTYPE inv_last1,
					  UDTYPE modulus_fac_last[3],
					  UDTYPE modulus_fac0[3],
					  UDTYPE modulus_fac1[3])
{

	UDTYPE last_temp;

	for (size_t j = 0; j < BRAMSIZE; j++)
	{
		for (size_t k = 0; k < 2; k++)
		{
			for (size_t i = 0; i < BRAMNUM; i += BNC)
			{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = encrypted0 inter false
#pragma HLS DEPENDENCE variable = encrypted1 inter false

				for (size_t c = 0; c < BNC; c++)
				{

					size_t i_c;
					i_c = i + c;

					last_temp = barrett_reduce_63(encrypted_last[k][i_c][j] + half_last, modulus_fac_last[0], modulus_fac_last[2]);

					rescale_core(encrypted0[k][i_c][j], last_temp, half_mod0, encrypted0[k][i_c][j],
								 modulus_fac0[0], modulus_fac0[1], modulus_fac0[2], inv_last0);
					rescale_core(encrypted1[k][i_c][j], last_temp, half_mod1, encrypted1[k][i_c][j],
								 modulus_fac1[0], modulus_fac1[1], modulus_fac1[2], inv_last1);
				}
			}
		}
	}
}

// BNC = BRAMNUM para count
template <unsigned BNC>
void rescale_para_mc1(UDTYPE encrypted_last[2][BRAMNUM][BRAMSIZE],
					  UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE],
					  UDTYPE half_last, UDTYPE half_mod0,
					  UDTYPE inv_last0,
					  UDTYPE modulus_fac_last[3],
					  UDTYPE modulus_fac0[3])
{

	UDTYPE last_temp;

	for (size_t j = 0; j < BRAMSIZE; j++)
	{
		for (size_t k = 0; k < 2; k++)
		{
			for (size_t i = 0; i < BRAMNUM; i += BNC)
			{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = encrypted0 inter false

				for (size_t c = 0; c < BNC; c++)
				{

					size_t i_c;
					i_c = i + c;

					last_temp = barrett_reduce_63(encrypted_last[k][i_c][j] + half_last, modulus_fac_last[0], modulus_fac_last[3]);

					rescale_core(encrypted0[k][i_c][j], last_temp, half_mod0, encrypted0[k][i_c][j],
								 modulus_fac0[0], modulus_fac0[1], modulus_fac0[2], inv_last0);
				}
			}
		}
	}
}

/*template <unsigned mod_count_in, unsigned mod_count_out>
void mod_switch_scale_to_next_new(UDTYPE encrypted[mod_count_in][2][BRAMNUM][BRAMSIZE],
								  UDTYPE rp[mod_count_in][N],
								  UDTYPE srp[mod_count_in][N],
								  UDTYPE idt[mod_count_in][N],
								  UDTYPE sidt[mod_count_in][N],
								  UDTYPE modulus[mod_count_in],
								  UDTYPE modulus_ratio[mod_count_in][2],
								  UDTYPE inv_last_coeff_mod_array[mod_count_out])

{
#pragma HLS INLINE off
	size_t coeff_mod_count = mod_count_in;
	size_t next_coeff_mod_count = mod_count_out;
	size_t encrypted_size = 2;

	UDTYPE encrypted_next_temp;
	UDTYPE half_mod[mod_count_out];
#pragma HLS ARRAY_PARTITION variable=half_mod complete dim=0

	UDTYPE last_modulus = modulus[next_coeff_mod_count];
	UDTYPE modulus_ratio1 = modulus_ratio[next_coeff_mod_count][1];
	UDTYPE half = last_modulus >> 1;

	for (int i = 0; i < next_coeff_mod_count; i++)
	{
		half_mod[i] = barrett_reduce_63(half, modulus[i], modulus_ratio[i][1]);
	}

intt_coeff_1_1:
	for (size_t j = 0; j < coeff_mod_count; j++)
	{
#pragma HLS UNROLL

#pragma HLS DEPENDENCE variable = encrypted inter false
		//#pragma HLS DEPENDENCE variable=idt inter false
		//#pragma HLS DEPENDENCE variable=sidt inter false

		intt_4core_mods<2, CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX>(
			encrypted[j], idt[j], sidt[j], modulus[j]);
	}


	for (size_t k = 0; k < encrypted_size; k++)
	{
	BRAMSIZE_loop:
		for (size_t j = 0; j < BRAMSIZE; j++)
		{
		BRAMNUM_loop:
			for (size_t i = 0; i < BRAMNUM; i++)
			{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = encrypted inter false

				encrypted_next_temp = barrett_reduce_63(encrypted[next_coeff_mod_count][k][i][j] + half, last_modulus, modulus_ratio1);

				for (size_t mod_index = 0; mod_index < next_coeff_mod_count; mod_index++)
				{
					rescale_core(encrypted[mod_index][k][i][j], encrypted_next_temp, half_mod[mod_index], encrypted[mod_index][k][i][j],
					modulus[mod_index], modulus_ratio[mod_index][0], modulus_ratio[mod_index][1], inv_last_coeff_mod_array[mod_index]);
				}
			}
		}
	}

	UDTYPE modulus_last = modulus[5];
	UDTYPE modulus_last_ratio[2];
	modulus_last_ratio[0] = modulus_ratio[5][0];
	modulus_last_ratio[1] = modulus_ratio[5][1];

	rescale_para_5<1>(encrypted[5], encrypted[0], encrypted[1], encrypted[2], encrypted[3], encrypted[4],
			half, half_mod, modulus_last, modulus_last_ratio, modulus, modulus_ratio, inv_last_coeff_mod_array);


	for (size_t j = 0; j < next_coeff_mod_count; j++)
	{
#pragma HLS UNROLL
		ntt_4core_mods<2, CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX>(
			encrypted[j], idt[j], sidt[j], modulus[j]);
	}
}*/

/*void rescale_mc3(UDTYPEin *in0,UDTYPEin *in1,UDTYPEin *in2,
		UDTYPE buf_last[2][BRAMNUM][BRAMSIZE], UDTYPE buf_m0[2][BRAMNUM][BRAMSIZE],
		UDTYPE buf_m1[2][BRAMNUM][BRAMSIZE], UDTYPE buf_m2[2][BRAMNUM][BRAMSIZE],
		UDTYPE mod_last[3], UDTYPE mod_m0[3], UDTYPE mod_m1[3], UDTYPE mod_m2[3],
		UDTYPE idt_last[N], UDTYPE idt_m0[N], UDTYPE idt_m1[N], UDTYPE idt_m2[N],
		UDTYPE sidt_last[N], UDTYPE sidt_m0[N], UDTYPE sidt_m1[N], UDTYPE sidt_m2[N],
		UDTYPE half_last, UDTYPE half_m0, UDTYPE half_m1, UDTYPE half_m2,
		UDTYPE inv_m0, UDTYPE inv_m1, UDTYPE inv_m2,
		UDTYPE rp_m0[N], UDTYPE rp_m1[N], UDTYPE rp_m2[N],
		UDTYPE srp_m0[N], UDTYPE srp_m1[N], UDTYPE srp_m2[N],
		UDTYPE out_m1[2][BRAMNUM][BRAMSIZE], UDTYPE out_m2[2][BRAMNUM][BRAMSIZE]){

	multiply_plain_in2_mc3(in0, in1, in2,
			buf_last, buf_m2, buf_m1,
			mod_last, mod_m2, mod_m1);

	in0 += N;
	in1 += N;
	in2 += N;

	intt_mc3<CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX>(
			buf_last, buf_m2, buf_m1,
			idt_last, idt_m2, idt_m1,
			sidt_last, sidt_m2, sidt_m1,
			mod_last[0], mod_m2[0], mod_m1[0]);

	rescale_para_mc3<1>(buf_last, buf_m2, buf_m1, buf_m0,
			half_last, half_m2, half_m1, half_m0,
			inv_m2, inv_m1, inv_m0,
			mod_last, mod_m2, mod_m1, mod_m0);

	ntt_mc3<CORENUM, BRAMNUM, BRAMSIZE, STAGEMAX>(
			buf_m2, buf_m1, buf_m0,
			rp_m2,  rp_m1, rp_m0,
			srp_m2,  srp_m1, srp_m0,
			mod_m2[0], mod_m1[0], mod_m0[0]);

	add_inplace_mc2<1>(out_m2, out_m1,
			buf_m2, buf_m1,
			mod_m2[0], mod_m1[0]);
}*/
