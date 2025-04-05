#include "define.h"

void multplain_inplace_cp_mc3(UDTYPEin *in0,UDTYPEin *in1,UDTYPEin *in2,
		UDTYPE result0[2][BRAMNUM][BRAMSIZE], UDTYPE result1[2][BRAMNUM][BRAMSIZE], UDTYPE result2[2][BRAMNUM][BRAMSIZE],
		UDTYPE modulus0[3], UDTYPE modulus1[3],UDTYPE modulus2[3]);

void multplain_p_mc2(UDTYPEin *in0,
		UDTYPE input0[2][BRAMNUM][BRAMSIZE], UDTYPE input1[2][BRAMNUM][BRAMSIZE],
		UDTYPE result0[2][BRAMNUM][BRAMSIZE], UDTYPE result1[2][BRAMNUM][BRAMSIZE],
		UDTYPE modulus0[3], UDTYPE modulus1[3]);

void multiply_plain_in1_mc2(UDTYPEin *in0,
		UDTYPE input0[2][BRAMNUM][BRAMSIZE], UDTYPE input1[2][BRAMNUM][BRAMSIZE],
		UDTYPE result0[2][BRAMNUM][BRAMSIZE], UDTYPE result1[2][BRAMNUM][BRAMSIZE],
		UDTYPE modulus0[3], UDTYPE modulus1[3]);

void multplain_inplace_cp_mc2(UDTYPEin *in0,UDTYPEin *in1,UDTYPEin *in_plain,
		UDTYPE result0[2][BRAMNUM][BRAMSIZE], UDTYPE result1[2][BRAMNUM][BRAMSIZE],
		UDTYPE modulus0[3], UDTYPE modulus1[3]);

void square_mc3(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE], UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE],
		UDTYPE encrypted2[2][BRAMNUM][BRAMSIZE], UDTYPE target0[BRAMNUM][BRAMSIZE],
		UDTYPE target1[BRAMNUM][BRAMSIZE], UDTYPE target2[BRAMNUM][BRAMSIZE],
	    UDTYPE modulus_fac0[3],UDTYPE modulus_fac1[3],UDTYPE modulus_fac2[3]);

/*void square_mc1(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE], UDTYPE target0[BRAMNUM][BRAMSIZE],
	    UDTYPE modulus_fac0[3]);*/


// ���ڲ��ж���1�ģ����Ե������ģ����С
template <unsigned modcount>
void square_mc1(UDTYPE encrypted[modcount][2][BRAMNUM][BRAMSIZE],
		UDTYPE target[modcount][BRAMNUM][BRAMSIZE],
	    UDTYPE modulus_fac[modcount][3])
{
#pragma HLS INLINE off

	UDTYPE operand0_0, operand0_1;
	UDTYPE operand1_0, operand1_1;
	UDTYPE operand2_0, operand2_1;

	UDTYPE result0_0, result0_1, result0_01;


	UDTYPE modulus0;
	UDTYPE ratio0;
	UDTYPE tratio0;

//	modulus0 = modulus_fac[0];
//	ratio0 = modulus_fac[1];
//	tratio0 = modulus_fac[2];


square_EE_loop:
	for (size_t m = 0; m < modcount; m++){

		modulus0 = modulus_fac[m][0];
		ratio0 = modulus_fac[m][1];
		tratio0 = modulus_fac[m][2];

		loop1:
		for(size_t j = 0; j < BRAMSIZE ;j++){
			loop2:
			for (size_t i = 0; i < BRAMNUM; i++){
//#pragma HLS UNROLL
#pragma HLS DEPENDENCE variable=target inter false
#pragma HLS DEPENDENCE variable=encrypted inter false
#pragma HLS PIPELINE

				//modulus0
				operand0_0 = encrypted[m][0][i][j];
				operand0_1 = encrypted[m][1][i][j];

				result0_0 =  HMult_single_return(operand0_0, operand0_0, modulus0, ratio0, tratio0);
				result0_1 =  HMult_single_return(operand0_1, operand0_1, modulus0, ratio0, tratio0);
				result0_01 =  HMult_single_return(operand0_0, operand0_1, modulus0, ratio0, tratio0);

				encrypted[m][0][i][j] = result0_0;
				encrypted[m][1][i][j] = result0_1;
				target[m][i][j] = Hadd_single(result0_01, result0_01, modulus0);

			}
		}
	}
}

template <unsigned paracount, unsigned modcount, unsigned bramnum, unsigned bramsize>
void pcmult2_inplace(
    UDTYPE encrypted1[modcount][2][bramnum][bramsize],
    UDTYPE plaintext[modcount][bramnum][bramsize],
    UDTYPE modulus0[5])
{
    #pragma HLS INLINE off


    add:
    for (uint32_t m = 0; m < modcount; m += paracount)
    {
        for (uint32_t j = 0; j < bramsize; j++)
        {
            for (uint32_t i = 0; i < bramnum; i++)
            {
                #pragma HLS DEPENDENCE variable=encrypted1 inter false
                #pragma HLS DEPENDENCE variable=plaintext inter false
                #pragma HLS PIPELINE

                for (uint32_t c = 0; c < paracount; c++)
                {
                    #pragma HLS UNROLL
                    uint32_t m_c = m + c;
                    UDTYPE temp = plaintext[m_c][i][j];
                    encrypted1[m_c][0][i][j] = HMult_single_return(encrypted1[m_c][0][i][j], temp, modulus0[0], modulus0[1], modulus0[2]);
                    encrypted1[m_c][1][i][j] = HMult_single_return(encrypted1[m_c][1][i][j], temp, modulus0[0], modulus0[1], modulus0[2]);
                }
            }
        }
    }
}


template <unsigned paracount, unsigned modcount, unsigned bramnum, unsigned bramsize>
void pcmult_inplace(
	UDTYPE encrypted1[modcount][2][bramnum][bramsize],
	UDTYPE plaintext[modcount][bramnum][bramsize],
	UDTYPE modulus0)
{
#pragma HLS INLINE off

	//UDTYPE modulus_temp;

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

					//modulus_temp = modulus[m_c];
					UDTYPE temp = plaintext[m_c][i][j];
					encrypted1[m_c][0][i][j] = HMult_single_return(encrypted1[m_c][0][i][j], temp, modulus0[0], modulus0[1], modulus0[2]);
					encrypted1[m_c][1][i][j] = HMult_single_return(encrypted1[m_c][1][i][j], temp, modulus0[0], modulus0[1], modulus0[2]);
				}
			}
		}
	}
}


/*template <unsigned mod_count>
void square(UDTYPE encrypted[mod_count][2][BRAMNUM][BRAMSIZE],
			   //UDTYPE destination[mod_count][2][bramnum][bramsize],
			   UDTYPE target[mod_count][BRAMNUM][BRAMSIZE],
			   UDTYPE modulus[mod_count],
			   UDTYPE modulus_ratio[mod_count][2])
{
#pragma HLS INLINE off

	UDTYPE operand0_0, operand0_1;
	UDTYPE operand1_0, operand1_1;
	UDTYPE operand2_0, operand2_1;

	UDTYPE result0_0, result0_1, result0_01;
	UDTYPE result1_0, result1_1, result1_01;
	UDTYPE result2_0, result2_1, result2_01;

	UDTYPE modulus0, modulus1, modulus2;
	UDTYPE ratio0, ratio1, ratio2;
	UDTYPE tratio0, tratio1, tratio2;

square_EE_loop:
	for (uint32_t m = 0; m < mod_count; m+=3)
	{
		for(size_t j = 0; j < BRAMSIZE ;j++){
			for (size_t i = 0; i < BRAMNUM; i++){
#pragma HLS DEPENDENCE variable=target inter false
#pragma HLS DEPENDENCE variable=encrypted inter false
#pragma HLS PIPELINE

				//modulus0
				modulus0 = modulus[m];
				ratio0 = modulus_ratio[m][0];
				tratio0 = modulus_ratio[m][1];
				operand0_0 = encrypted[m][0][i][j];
				operand0_1 = encrypted[m][1][i][j];

				result0_0 =  HMult_single_return(operand0_0, operand0_0, modulus0, ratio0, tratio0);
				result0_1 =  HMult_single_return(operand0_1, operand0_1, modulus0, ratio0, tratio0);
				result0_01 =  HMult_single_return(operand0_0, operand0_1, modulus0, ratio0, tratio0);

				encrypted[m][0][i][j] = result0_0;
				encrypted[m][1][i][j] = result0_1;
				target[m][i][j] = Hadd_single(result0_01, result0_01, modulus0);

				//modulus1
				modulus1 = modulus[m + 1];
				ratio1 = modulus_ratio[m + 1][0];
				tratio1 = modulus_ratio[m + 1][1];
				operand1_0 = encrypted[m + 1][0][i][j];
				operand1_1 = encrypted[m + 1][1][i][j];

				result1_0 =  HMult_single_return(operand1_0, operand1_0, modulus1, ratio1, tratio1);
				result1_1 =  HMult_single_return(operand1_1, operand1_1, modulus1, ratio1, tratio1);
				result1_01 =  HMult_single_return(operand1_0, operand1_1, modulus1, ratio1, tratio1);

				encrypted[m + 1][0][i][j] = result1_0;
				encrypted[m + 1][1][i][j] = result1_1;
				target[m + 1][i][j] = Hadd_single(result1_01, result1_01, modulus1);

				//modulus2
				modulus2 = modulus[m + 2];
				ratio2 = modulus_ratio[m + 2][0];
				tratio2 = modulus_ratio[m + 2][1];
				operand2_0 = encrypted[m + 2][0][i][j];
				operand2_1 = encrypted[m + 2][1][i][j];

				result2_0 =  HMult_single_return(operand2_0, operand2_0, modulus2, ratio2, tratio2);
				result2_1 =  HMult_single_return(operand2_1, operand2_1, modulus2, ratio2, tratio2);
				result2_01 =  HMult_single_return(operand2_0, operand2_1, modulus2, ratio2, tratio2);

				encrypted[m + 2][0][i][j] = result2_0;
				encrypted[m + 2][1][i][j] = result2_1;
				target[m + 2][i][j] = Hadd_single(result2_01, result2_01, modulus2);
			}
		}
	}
}*/

/*template <unsigned bramnum, unsigned bramsize>
void multiply_plain_axi2_single(UDTYPEin *in, UDTYPE result[2][bramnum][bramsize],
							UDTYPE modulus, UDTYPE ratio_0, UDTYPE ratio_1)
{
#pragma HLS INLINE off

	UDTYPEin in_temp;
	UDTYPE plain, encrypted_0, encrypted_1;

Hmult_plain_j:
	for (uint32_t j = 0; j < bramsize; j++)
	{
		for (uint32_t i = 0; i < bramnum; i++)
		{
#pragma HLS PIPELINE
			in_temp = *in;
			in++;

			encrypted_0 = (in_temp)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
			encrypted_1 = (in_temp)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
			plain = (in_temp)((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);

			result[0][i][j] = HMult_single_return(encrypted_0, plain, modulus, ratio_0, ratio_1);
			result[1][i][j] = HMult_single_return(encrypted_1, plain, modulus, ratio_0, ratio_1);
		}
	}
}*/

/*template <unsigned mod_count, unsigned bramnum, unsigned bramsize>
void multiply_plain_axi2(UDTYPEin *in0,
				 UDTYPE result[mod_count][2][bramnum][bramsize],
				 UDTYPE modulus[INITMODCOUNT],
				 UDTYPE modulus_ratio[INITMODCOUNT][2])

{
#pragma HLS INLINE off

Hmult_plain_m:
	for (uint32_t m = 0; m < 6; m++)
	{
		multiply_plain_axi2_single<bramnum, bramsize>(in0, result[m], modulus[m], modulus_ratio[m][0], modulus_ratio[m][1]);
		in0 += N;
	}
}*/

/*template <unsigned mod_count, unsigned bramnum, unsigned bramsize>
void multiply_plain_axi(
		UDTYPEin *in0, UDTYPEin *in1, UDTYPEin *in2,
		UDTYPE encrypted[mod_count][2][bramnum][bramsize],
		UDTYPE modulus[INITMODCOUNT],
		UDTYPE modulus_ratio[INITMODCOUNT][2])
{
#pragma HLS INLINE off

	UDTYPEin in_temp0, in_temp1, in_temp2;
	UDTYPE plain0, plain1, plain2;
	UDTYPE encrypted0_0, encrypted0_1;
	UDTYPE encrypted1_0, encrypted1_1;
	UDTYPE encrypted2_0, encrypted2_1;
	UDTYPE modulus0, modulus1, modulus2;
	UDTYPE ratio0, ratio1, ratio2;
	UDTYPE tratio0, tratio1, tratio2;

	Hmult_plain_m:
	for (uint32_t m = 0; m < mod_count; m+=3)
	{
		Hmult_plain_j:
		for (uint32_t j = 0; j < bramsize; j++)
		{
			for (uint32_t i = 0; i < bramnum; i++)
			{
#pragma HLS DEPENDENCE variable=encrypted inter false
#pragma HLS PIPELINE
				in_temp0 = *in0;
				in0++;

				modulus0 = modulus[m];
				ratio0 = modulus_ratio[m][0];
				tratio0 = modulus_ratio[m][1];
				plain0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
				encrypted0_0 = encrypted[m][0][i][j];
				encrypted0_1 = encrypted[m][1][i][j];

				encrypted[m][0][i][j] = HMult_single_return(encrypted0_0, plain0, modulus0, ratio0, tratio0);
				encrypted[m][1][i][j] = HMult_single_return(encrypted0_1, plain0, modulus0, ratio0, tratio0);

				in_temp1 = *in1;
				in1++;

				modulus1 = modulus[m + 1];
				ratio1 = modulus_ratio[m + 1][0];
				tratio1 = modulus_ratio[m + 1][1];
				plain1 = (in_temp1)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
				encrypted1_0 = encrypted[m + 1][0][i][j];
				encrypted1_1 = encrypted[m + 1][1][i][j];

				encrypted[m + 1][0][i][j] = HMult_single_return(encrypted1_0, plain1, modulus1, ratio1, tratio1);
				encrypted[m + 1][1][i][j] = HMult_single_return(encrypted1_1, plain1, modulus1, ratio1, tratio1);

				in_temp2 = *in2;
				in2++;

				modulus2 = modulus[m + 2];
				ratio2 = modulus_ratio[m + 2][0];
				tratio2 = modulus_ratio[m + 2][1];
				plain2 = (in_temp2)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
				encrypted2_0 = encrypted[m + 2][0][i][j];
				encrypted2_1 = encrypted[m + 2][1][i][j];

				encrypted[m + 2][0][i][j] = HMult_single_return(encrypted2_0, plain2, modulus2, ratio2, tratio2);
				encrypted[m + 2][1][i][j] = HMult_single_return(encrypted2_1, plain2, modulus2, ratio2, tratio2);
			}
		}
	}
}*/


/*template <unsigned mod_count, unsigned bramnum, unsigned bramsize>
void multiply_plain_axi_port1(
		UDTYPEin *in0,
		UDTYPE in[mod_count][2][bramnum][bramsize],
		UDTYPE out[mod_count][2][bramnum][bramsize],
		UDTYPE modulus[INITMODCOUNT],
		UDTYPE modulus_ratio[INITMODCOUNT][2])
{
#pragma HLS INLINE off

	UDTYPEin in_temp0, in_temp1, in_temp2;
	UDTYPE plain0, plain1, plain2;
	UDTYPE encrypted0_0, encrypted0_1;
	UDTYPE encrypted1_0, encrypted1_1;
	UDTYPE encrypted2_0, encrypted2_1;
	UDTYPE modulus0, modulus1, modulus2;
	UDTYPE ratio0, ratio1, ratio2;
	UDTYPE tratio0, tratio1, tratio2;

	Hmult_plain_m:
	for (uint32_t m = 0; m < mod_count; m+=2)
	{
		Hmult_plain_j:
		for (uint32_t j = 0; j < bramsize; j++)
		{
			for (uint32_t i = 0; i < bramnum; i++)
			{
#pragma HLS DEPENDENCE variable=in inter false
#pragma HLS DEPENDENCE variable=out inter false
#pragma HLS PIPELINE
				in_temp0 = *in0;
				in0++;

				modulus0 = modulus[m];
				ratio0 = modulus_ratio[m][0];
				tratio0 = modulus_ratio[m][1];
				plain0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
				encrypted0_0 = in[m][0][i][j];
				encrypted0_1 = in[m][1][i][j];

				out[m][0][i][j] = HMult_single_return(encrypted0_0, plain0, modulus0, ratio0, tratio0);
				out[m][1][i][j] = HMult_single_return(encrypted0_1, plain0, modulus0, ratio0, tratio0);

				modulus1 = modulus[m + 1];
				ratio1 = modulus_ratio[m + 1][0];
				tratio1 = modulus_ratio[m + 1][1];
				plain1 = (in_temp0)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
				encrypted1_0 = in[m + 1][0][i][j];
				encrypted1_1 = in[m + 1][1][i][j];

				out[m + 1][0][i][j] = HMult_single_return(encrypted1_0, plain1, modulus1, ratio1, tratio1);
				out[m + 1][1][i][j] = HMult_single_return(encrypted1_1, plain1, modulus1, ratio1, tratio1);

			}
		}
	}
}*/







