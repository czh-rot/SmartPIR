
#include "define.h"

/*template <unsigned polys_size, unsigned mod_count>
void add_inplace(
	UDTYPE encrypted1[2][mod_count][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted2[2][mod_count][BRAMNUM][BRAMSIZE],
	UDTYPE coeff_modulus[ROOTMODCOUNT] // encrypted1 + encrypted2�ѽ���浽encrypted1
)
{

	size_t coeff_count = N;

	for (size_t j = 0; j < 2; j++)
	{
		for (size_t i = 0; i < mod_count; i++)
		{
			add_8_1024(encrypted1[j][i],
					   encrypted2[j][i], coeff_modulus[i],
					   encrypted1[j][i]);
		}
	}
}*/

template <unsigned paracount, unsigned modcount, unsigned bramnum, unsigned bramsize>
void add_inplace(
	UDTYPE encrypted1[modcount][2][bramnum][bramsize],
	UDTYPE encrypted2[modcount][2][bramnum][bramsize],
	UDTYPE modulus[5])
{
#pragma HLS INLINE off

	UDTYPE modulus_temp;

    add:
	for (uint32_t m = 0; m < modcount; m += paracount)
	{
		for (uint32_t j = 0; j < bramsize; j++)
		{
			for (uint32_t i = 0; i < bramnum; i++)
			{

#pragma HLS DEPENDENCE variable = encrypted1 inter false
#pragma HLS DEPENDENCE variable = encrypted2 inter false
#pragma HLS PIPELINE

				for (uint32_t c = 0; c < paracount; c++)
				{
					uint32_t m_c;
					m_c = m + c;

					modulus_temp = modulus[m_c];

					encrypted1[m_c][0][i][j] = Hadd_single(encrypted1[m_c][0][i][j], encrypted2[m_c][0][i][j], modulus_temp);
					encrypted1[m_c][1][i][j] = Hadd_single(encrypted1[m_c][1][i][j], encrypted2[m_c][1][i][j], modulus_temp);
				}
			}
		}
	}
}

template <unsigned paracount, unsigned modcount, unsigned bramnum, unsigned bramsize>
void add_plain_inplace(
	UDTYPE encrypted1[modcount][2][bramnum][bramsize],
	UDTYPE plaintext[modcount][bramnum][bramsize],
	UDTYPE modulus[INITMODCOUNT])
{
#pragma HLS INLINE off

	UDTYPE modulus_temp;

    add:
	for (uint32_t m = 0; m < modcount; m += paracount)
	{
		for (uint32_t j = 0; j < bramsize; j++)
		{
			for (uint32_t i = 0; i < bramnum; i++)
			{

#pragma HLS DEPENDENCE variable = encrypted1 inter false
#pragma HLS DEPENDENCE variable = plaintext inter false
#pragma HLS PIPELINE

				for (uint32_t c = 0; c < paracount; c++)
				{
					uint32_t m_c;
					m_c = m + c;

					modulus_temp = modulus[m_c];

					UDTYPE tmp = plaintext[m_c][i][j];
					encrypted1[m_c][0][i][j] = Hadd_single(encrypted1[m_c][0][i][j], tmp, modulus_temp);
					encrypted1[m_c][1][i][j] = Hadd_single(encrypted1[m_c][1][i][j], tmp, modulus_temp);
				}
			}
		}
	}
}


//BNC = BRAMNUM para count
template <unsigned BNC>
void add_inplace_mc3(
	UDTYPE encrypted_save0[2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted_save1[2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted_save2[2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted_in0[2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted_in1[2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted_in2[2][BRAMNUM][BRAMSIZE],
	UDTYPE modulus0,UDTYPE modulus1,UDTYPE modulus2)
{
#pragma HLS INLINE off

	UDTYPE modulus_temp;

    add:

	for (uint32_t j = 0; j < BRAMSIZE; j++)
	{
		for (uint32_t k = 0; k < 2; k++)
		{
			for (uint32_t i = 0; i < BRAMNUM; i += BNC)
			{

#pragma HLS DEPENDENCE variable = encrypted_save0 inter false
#pragma HLS DEPENDENCE variable = encrypted_in0 inter false
#pragma HLS DEPENDENCE variable = encrypted_save1 inter false
#pragma HLS DEPENDENCE variable = encrypted_in1 inter false
#pragma HLS DEPENDENCE variable = encrypted_save2 inter false
#pragma HLS DEPENDENCE variable = encrypted_in2 inter false
#pragma HLS PIPELINE

				for (uint32_t c = 0; c < BNC; c++)
				{
					uint32_t i_c;
					i_c = i + c;

					encrypted_save0[k][i_c][j] = Hadd_single(encrypted_save0[k][i_c][j], encrypted_in0[k][i_c][j], modulus0);
					encrypted_save1[k][i_c][j] = Hadd_single(encrypted_save1[k][i_c][j], encrypted_in1[k][i_c][j], modulus1);
					encrypted_save2[k][i_c][j] = Hadd_single(encrypted_save2[k][i_c][j], encrypted_in2[k][i_c][j], modulus2);

				}
			}
		}
	}
}

//BNC = BRAMNUM para count
template <unsigned BNC>
void add_inplace_mc2(
	UDTYPE encrypted_save0[2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted_save1[2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted_in0[2][BRAMNUM][BRAMSIZE],
	UDTYPE encrypted_in1[2][BRAMNUM][BRAMSIZE],
	UDTYPE modulus0, UDTYPE modulus1)
{
#pragma HLS INLINE off

	UDTYPE modulus_temp;

    add:

	for (uint32_t j = 0; j < BRAMSIZE; j++)
	{
		for (uint32_t k = 0; k < 2; k++)
		{
			for (uint32_t i = 0; i < BRAMNUM; i += BNC)
			{

#pragma HLS DEPENDENCE variable = encrypted_save0 inter false
#pragma HLS DEPENDENCE variable = encrypted_in0 inter false
#pragma HLS DEPENDENCE variable = encrypted_save1 inter false
#pragma HLS DEPENDENCE variable = encrypted_in1 inter false
#pragma HLS PIPELINE

				for (uint32_t c = 0; c < BNC; c++)
				{
					uint32_t i_c;
					i_c = i + c;

					encrypted_save0[k][i_c][j] = Hadd_single(encrypted_save0[k][i_c][j], encrypted_in0[k][i_c][j], modulus0);
					encrypted_save1[k][i_c][j] = Hadd_single(encrypted_save1[k][i_c][j], encrypted_in1[k][i_c][j], modulus1);

				}
			}
		}
	}
}
