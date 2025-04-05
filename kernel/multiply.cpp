#include "top.h"

using namespace std;

void multplain_inplace_cp_mc3(UDTYPEin *in0, UDTYPEin *in1, UDTYPEin *in2,
							  UDTYPE result0[2][BRAMNUM][BRAMSIZE], UDTYPE result1[2][BRAMNUM][BRAMSIZE], UDTYPE result2[2][BRAMNUM][BRAMSIZE],
							  UDTYPE modulus0[3], UDTYPE modulus1[3], UDTYPE modulus2[3])
{
#pragma HLS INLINE off

	UDTYPEin in_temp[3];
	UDTYPEin encrypted[3][2];
	UDTYPEin plain[3];
#pragma HLS ARRAY_PARTITION variable = in_temp complete dim = 0
#pragma HLS ARRAY_PARTITION variable = encrypted complete dim = 0
#pragma HLS ARRAY_PARTITION variable = plain complete dim = 0

	IDXTYPE i, j;

Hmult_loop:
	for (IDXTYPE idx = 0; idx < N; idx++)
	{
#pragma HLS PIPELINE
		// i = idx % BRAMNUM;
		// j = (idx / BRAMNUM) % BRAMSIZE;

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);

		in_temp[0] = in0[idx];
		in_temp[1] = in1[idx];
		in_temp[2] = in2[idx];

		encrypted[0][0] = (in_temp[0])((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		encrypted[0][1] = (in_temp[0])((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		plain[0] = (in_temp[0])((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);

		// test
		// cout << encrypted[0][0] << endl;
		// cout << encrypted[0][1] << endl;
		// cout << plain[0] << endl;

		result0[0][i][j] = HMult_single_return(encrypted[0][0], plain[0], modulus0[0], modulus0[1], modulus0[2]);
		result0[1][i][j] = HMult_single_return(encrypted[0][1], plain[0], modulus0[0], modulus0[1], modulus0[2]);

		// cout << "result0[0][i][j] = " << result0[0][i][j] << endl;
		// cout << "result0[1][i][j] = " << result0[1][i][j] << endl;

		encrypted[1][0] = (in_temp[1])((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		encrypted[1][1] = (in_temp[1])((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		plain[1] = (in_temp[1])((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);
		result1[0][i][j] = HMult_single_return(encrypted[1][0], plain[1], modulus1[0], modulus1[1], modulus1[2]);
		result1[1][i][j] = HMult_single_return(encrypted[1][1], plain[1], modulus1[0], modulus1[1], modulus1[2]);

		encrypted[2][0] = (in_temp[2])((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		encrypted[2][1] = (in_temp[2])((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		plain[2] = (in_temp[2])((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);
		result2[0][i][j] = HMult_single_return(encrypted[2][0], plain[2], modulus2[0], modulus2[1], modulus2[2]);
		result2[1][i][j] = HMult_single_return(encrypted[2][1], plain[2], modulus2[0], modulus2[1], modulus2[2]);
	}
}

void multplain_inplace_cp_mc2(UDTYPEin *in0, UDTYPEin *in1, UDTYPEin *in_plain,
							  UDTYPE result0[2][BRAMNUM][BRAMSIZE], UDTYPE result1[2][BRAMNUM][BRAMSIZE],
							  UDTYPE modulus0[3], UDTYPE modulus1[3])
{
#pragma HLS INLINE off

	UDTYPEin in_temp[3];
	UDTYPEin encrypted[3][2];
	UDTYPEin plain[3];
#pragma HLS ARRAY_PARTITION variable = in_temp complete dim = 0
#pragma HLS ARRAY_PARTITION variable = encrypted complete dim = 0
#pragma HLS ARRAY_PARTITION variable = plain complete dim = 0

	IDXTYPE i, j;

Hmult_loop:
	for (IDXTYPE idx = 0; idx < N; idx++)
	{
#pragma HLS PIPELINE
		// i = idx % BRAMNUM;
		// j = (idx / BRAMNUM) % BRAMSIZE;

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);

		in_temp[0] = in0[idx];
		in_temp[1] = in1[idx];
		in_temp[2] = in_plain[idx];

		plain[0] = (in_temp[2])((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		plain[1] = (in_temp[2])((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);

		encrypted[0][0] = (in_temp[0])((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		encrypted[0][1] = (in_temp[0])((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		result0[0][i][j] = HMult_single_return(encrypted[0][0], plain[0], modulus0[0], modulus0[1], modulus0[2]);
		result0[1][i][j] = HMult_single_return(encrypted[0][1], plain[0], modulus0[0], modulus0[1], modulus0[2]);

		encrypted[1][0] = (in_temp[1])((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		encrypted[1][1] = (in_temp[1])((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		result1[0][i][j] = HMult_single_return(encrypted[1][0], plain[1], modulus1[0], modulus1[1], modulus1[2]);
		result1[1][i][j] = HMult_single_return(encrypted[1][1], plain[1], modulus1[0], modulus1[1], modulus1[2]);
	}
}

// �����in0����ʾ������ģ������
void multplain_p_mc2(UDTYPEin *in0,
					 UDTYPE input0[2][BRAMNUM][BRAMSIZE], UDTYPE input1[2][BRAMNUM][BRAMSIZE],
					 UDTYPE result0[2][BRAMNUM][BRAMSIZE], UDTYPE result1[2][BRAMNUM][BRAMSIZE],
					 UDTYPE modulus0[3], UDTYPE modulus1[3])
{
#pragma HLS INLINE off

	UDTYPEin in_temp; // ��Ϊ��3��ģ����������3�������������Ҫ�ֶ���
	UDTYPEin encrypted[2][2];
	UDTYPEin plain[2];
#pragma HLS ARRAY_PARTITION variable = encrypted complete dim = 0
#pragma HLS ARRAY_PARTITION variable = plain complete dim = 0

	IDXTYPE i, j;

Hmult_plain_j:
	for (IDXTYPE idx = 0; idx < N; idx++)
	{
#pragma HLS PIPELINE

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);

		in_temp = in0[idx];
		plain[0] = (in_temp)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		plain[1] = (in_temp)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);

		encrypted[0][0] = input0[0][i][j];
		encrypted[0][1] = input0[1][i][j];
		result0[0][i][j] = HMult_single_return(encrypted[0][0], plain[0], modulus0[0], modulus0[1], modulus0[2]);
		result0[1][i][j] = HMult_single_return(encrypted[0][1], plain[0], modulus0[0], modulus0[1], modulus0[2]);

		encrypted[1][0] = input1[0][i][j];
		encrypted[1][1] = input1[1][i][j];
		result1[0][i][j] = HMult_single_return(encrypted[1][0], plain[1], modulus1[0], modulus1[1], modulus1[2]);
		result1[1][i][j] = HMult_single_return(encrypted[1][1], plain[1], modulus1[0], modulus1[1], modulus1[2]);
	}
}

void square_mc3(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE], UDTYPE encrypted1[2][BRAMNUM][BRAMSIZE],
				UDTYPE encrypted2[2][BRAMNUM][BRAMSIZE], UDTYPE target0[BRAMNUM][BRAMSIZE],
				UDTYPE target1[BRAMNUM][BRAMSIZE], UDTYPE target2[BRAMNUM][BRAMSIZE],
				UDTYPE modulus_fac0[3], UDTYPE modulus_fac1[3], UDTYPE modulus_fac2[3])
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

	modulus0 = modulus_fac0[0];
	ratio0 = modulus_fac0[1];
	tratio0 = modulus_fac0[2];

	modulus1 = modulus_fac1[0];
	ratio1 = modulus_fac1[1];
	tratio1 = modulus_fac1[2];

	modulus2 = modulus_fac2[0];
	ratio2 = modulus_fac2[1];
	tratio2 = modulus_fac2[2];

square_EE_loop:

	for (size_t j = 0; j < BRAMSIZE; j++)
	{
		for (size_t i = 0; i < BRAMNUM; i++)
		{
#pragma HLS DEPENDENCE variable = target0 inter false
#pragma HLS DEPENDENCE variable = encrypted0 inter false
#pragma HLS DEPENDENCE variable = target1 inter false
#pragma HLS DEPENDENCE variable = encrypted1 inter false
#pragma HLS DEPENDENCE variable = target2 inter false
#pragma HLS DEPENDENCE variable = encrypted2 inter false
#pragma HLS PIPELINE

			// modulus0
			operand0_0 = encrypted0[0][i][j];
			operand0_1 = encrypted0[1][i][j];

			result0_0 = HMult_single_return(operand0_0, operand0_0, modulus0, ratio0, tratio0);
			result0_1 = HMult_single_return(operand0_1, operand0_1, modulus0, ratio0, tratio0);
			result0_01 = HMult_single_return(operand0_0, operand0_1, modulus0, ratio0, tratio0);

			encrypted0[0][i][j] = result0_0;
			encrypted0[1][i][j] = result0_1;
			target0[i][j] = Hadd_single(result0_01, result0_01, modulus0);

			// modulus1
			operand1_0 = encrypted1[0][i][j];
			operand1_1 = encrypted1[1][i][j];

			result1_0 = HMult_single_return(operand1_0, operand1_0, modulus1, ratio1, tratio1);
			result1_1 = HMult_single_return(operand1_1, operand1_1, modulus1, ratio1, tratio1);
			result1_01 = HMult_single_return(operand1_0, operand1_1, modulus1, ratio1, tratio1);

			encrypted1[0][i][j] = result1_0;
			encrypted1[1][i][j] = result1_1;
			target1[i][j] = Hadd_single(result1_01, result1_01, modulus1);

			// modulus2
			operand2_0 = encrypted2[0][i][j];
			operand2_1 = encrypted2[1][i][j];

			result2_0 = HMult_single_return(operand2_0, operand2_0, modulus2, ratio2, tratio2);
			result2_1 = HMult_single_return(operand2_1, operand2_1, modulus2, ratio2, tratio2);
			result2_01 = HMult_single_return(operand2_0, operand2_1, modulus2, ratio2, tratio2);

			encrypted2[0][i][j] = result2_0;
			encrypted2[1][i][j] = result2_1;
			target2[i][j] = Hadd_single(result2_01, result2_01, modulus2);
		}
	}
}

void multiply_plain_in1_mc2(UDTYPEin *in0,
							UDTYPE input0[2][BRAMNUM][BRAMSIZE], UDTYPE input1[2][BRAMNUM][BRAMSIZE],
							UDTYPE result0[2][BRAMNUM][BRAMSIZE], UDTYPE result1[2][BRAMNUM][BRAMSIZE],
							UDTYPE modulus0[3], UDTYPE modulus1[3])
{
#pragma HLS INLINE off

	UDTYPEin in_temp0, in_temp1, in_temp2;
	UDTYPE plain0, plain1, plain2;
	UDTYPE encrypted0_0, encrypted0_1;
	UDTYPE encrypted1_0, encrypted1_1;
	UDTYPE encrypted2_0, encrypted2_1;

	IDXTYPE i, j;

Hmult_plain_j:
	for (IDXTYPE idx = 0; idx < N; idx++)
	{
#pragma HLS PIPELINE

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);

		in_temp0 = in0[idx];
		plain0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		plain1 = (in_temp0)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);

		encrypted0_0 = input0[0][i][j];
		encrypted0_1 = input0[1][i][j];
		result0[0][i][j] = HMult_single_return(encrypted0_0, plain0, modulus0[0], modulus0[1], modulus0[2]);
		result0[1][i][j] = HMult_single_return(encrypted0_1, plain0, modulus0[0], modulus0[1], modulus0[2]);

		encrypted1_0 = input1[0][i][j];
		encrypted1_1 = input1[1][i][j];
		result1[0][i][j] = HMult_single_return(encrypted1_0, plain1, modulus1[0], modulus1[1], modulus1[2]);
		result1[1][i][j] = HMult_single_return(encrypted1_1, plain1, modulus1[0], modulus1[1], modulus1[2]);
	}
}

/*void square_mc1(UDTYPE encrypted0[2][BRAMNUM][BRAMSIZE], UDTYPE target0[BRAMNUM][BRAMSIZE],
		UDTYPE modulus_fac0[3])
{
#pragma HLS INLINE off

	UDTYPE operand0_0, operand0_1;
	UDTYPE operand1_0, operand1_1;
	UDTYPE operand2_0, operand2_1;

	UDTYPE result0_0, result0_1, result0_01;


	UDTYPE modulus0;
	UDTYPE ratio0;
	UDTYPE tratio0;

	modulus0 = modulus_fac0[0];
	ratio0 = modulus_fac0[1];
	tratio0 = modulus_fac0[2];


square_EE_loop:

	for(size_t j = 0; j < BRAMSIZE ;j++){
		for (size_t i = 0; i < BRAMNUM; i++){
#pragma HLS DEPENDENCE variable=target0 inter false
#pragma HLS DEPENDENCE variable=encrypted0 inter false
#pragma HLS PIPELINE

			//modulus0
			operand0_0 = encrypted0[0][i][j];
			operand0_1 = encrypted0[1][i][j];

			result0_0 =  HMult_single_return(operand0_0, operand0_0, modulus0, ratio0, tratio0);
			result0_1 =  HMult_single_return(operand0_1, operand0_1, modulus0, ratio0, tratio0);
			result0_01 =  HMult_single_return(operand0_0, operand0_1, modulus0, ratio0, tratio0);

			encrypted0[0][i][j] = result0_0;
			encrypted0[1][i][j] = result0_1;
			target0[i][j] = Hadd_single(result0_01, result0_01, modulus0);

		}
	}
}*/

/*void Hmult_N_ddr(UDTYPEin* operand1,
	UDTYPEin* operand2,
	UDTYPE modulus, UDTYPE const_ratio[2],UDTYPE result[BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	UDTYPE op1[PACKNUM_PMUL];
	UDTYPE op2[PACKNUM_PMUL];

	size_t index=0;
	size_t bramindex=0;

	UDTYPEin temp1,temp2;

	ddr_loop_j:
	for(size_t j = 0; j < BRAMSIZE ;j++){
		ddr_loop_i:
		for (size_t n = 0; n < BRAMNUM/(PACKNUM_PMUL); n++){
#pragma HLS PIPELINE
			temp1=*operand1;
			temp2=*operand2;
			operand1++;
			operand2++;

				   for(int m=0;m<PACKNUM_PMUL;m++){


					op1[m]=(temp1)((m+1)*BITWIDTH-1,m*BITWIDTH);
					op2[m]=(temp2)((m+1)*BITWIDTH-1,m*BITWIDTH);


				   result[bramindex][j] = HMult_single(op1[m],op2[m],modulus,const_ratio);
				   bramindex++;

			}
		}
	}
}*/

/*void Hmult_N_ddr(UDTYPEin* encrypted_in,
	UDTYPEin* plain_in,
	UDTYPE modulus,
	UDTYPE modulus_ratio_0, UDTYPE modulus_ratio_1,
	UDTYPE result[2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	UDTYPE encrypted[2][2];
#pragma HLS ARRAY_PARTITION variable=encrypted complete dim=1
#pragma HLS ARRAY_PARTITION variable=encrypted complete dim=2
	UDTYPE plain[2];
#pragma HLS ARRAY_PARTITION variable=plain complete dim=1

	UDTYPEin encrypted_temp, plain_temp;

	ddr_loop_j:
	for(size_t j = 0; j < BRAMSIZE ;j++){
		ddr_loop_i:
		for (size_t n = 0; n < BRAMNUM; n += 2){
#pragma HLS PIPELINE
			encrypted_temp = *encrypted_in;
			plain_temp = *plain_in;
			encrypted_in++;
			plain_in++;

			plain[0] = (plain_temp) ( (0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
			plain[1] = (plain_temp) ( (1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);

			encrypted[0][0] = (encrypted_temp) ( (0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
			encrypted[1][0] = (encrypted_temp) ( (1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
			encrypted[0][1] = (encrypted_temp) ( (2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);
			encrypted[1][1] = (encrypted_temp) ( (3 + 1) * BITWIDTH - 1, 3 * BITWIDTH);

			result[0][n][j] =  HMult_single_return(encrypted[0][0], plain[0], modulus, modulus_ratio_0, modulus_ratio_1);
			result[1][n][j] =  HMult_single_return(encrypted[1][0], plain[0], modulus, modulus_ratio_0, modulus_ratio_1);
			result[0][n + 1][j] =  HMult_single_return(encrypted[0][1], plain[1], modulus, modulus_ratio_0, modulus_ratio_1);
			result[1][n + 1][j] =  HMult_single_return(encrypted[1][1], plain[1], modulus, modulus_ratio_0, modulus_ratio_1);
		}
	}
}*/

/*void dyadic_product_coeffmod(UDTYPE operand1[BRAMNUM][BRAMSIZE],
	UDTYPE operand2[BRAMSIZE][BRAMSIZE],
	UDTYPE modulus, UDTYPE const_ratio[2], UDTYPE result[BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	coff_loop:
	for(size_t j = 0; j < BRAMSIZE ;j++){
#pragma HLS PIPELINE
		for (size_t i = 0; i < BRAMNUM; i++){
			UDTYPE z[2], tmp1, tmp2[2], tmp3, carry;
			UDTYPE2 tmpz, tmp4;
			tmpz = operand1[i][j] * operand2[i][j];
			z[0] = tmpz(BITWIDTH-1,0);
			z[1] = tmpz(2*BITWIDTH-1,BITWIDTH);

			// Round 1
			carry = (z[0] * const_ratio[0])(2*BITWIDTH-1,BITWIDTH);
			tmp4 = z[0] * const_ratio[1];
			tmp2[0] = tmp4(BITWIDTH-1,0);
			tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
			tmp3 = tmp2[1] + add_uint64(tmp2[0], carry, &tmp1);

			// Round 2
			tmp4 = z[1] * const_ratio[0];
			tmp2[0] = tmp4(BITWIDTH-1,0);
			tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
			carry = tmp2[1] + add_uint64(tmp1, tmp2[0], &tmp1);

			// This is all we care about
			tmp1 = z[1] * const_ratio[1] + tmp3 + carry;

			// Barrett subtraction
			tmp3 = z[0] - tmp1 * modulus;

			// Claim: One more subtraction is enough
			result[i][j] = tmp3 - (modulus & UDTYPE(
				-DTYPE(tmp3 >= modulus)));
		}
	}
}*/

/*void dyadic_product_coeffmod_cross(UDTYPE operand1[BRAMNUM][BRAMSIZE],
	UDTYPE operand2[BRAMSIZE][BRAMSIZE],
	UDTYPE operand1_[BRAMNUM][BRAMSIZE],
	UDTYPE operand2_[BRAMSIZE][BRAMSIZE],
	UDTYPE modulus, UDTYPE const_ratio[2], UDTYPE result[BRAMNUM][BRAMSIZE])
{

	UDTYPE result1,result2;
	coff_loop:for(size_t j = 0; j < BRAMSIZE ;j++){

		for (size_t i = 0; i < BRAMNUM; i++){
#pragma HLS PIPELINE
			HMult_single(operand1[i][j],operand2[i][j],result1,modulus,const_ratio); // @suppress("Invalid arguments")
			HMult_single(operand1_[i][j],operand2_[i][j],result2,modulus,const_ratio); // @suppress("Invalid arguments")

			result[i][j]=result1+result2;


		}
	}
}*/

/*void square_N(UDTYPE operand[BRAMNUM][BRAMSIZE],
	UDTYPE result[BRAMSIZE][BRAMSIZE],
	UDTYPE modulus, UDTYPE const_ratio[2])
{
#pragma HLS INLINE off

	UDTYPE operand_temp;
	UDTYPE result_temp;

	square_N_j:
	for(size_t j = 0; j < BRAMSIZE ;j++){
		square_N_i:
		for (size_t i = 0; i < BRAMNUM; i++){
#pragma HLS PIPELINE II=1
			operand_temp = operand[i][j];
			result_temp =  HMult_single_return(operand_temp, operand_temp, modulus, const_ratio);
			result[i][j] = result_temp;
		}
	}
}*/

/*void square_N_cross(UDTYPE operand1[BRAMNUM][BRAMSIZE], UDTYPE operand2[BRAMNUM][BRAMSIZE],
	UDTYPE result[BRAMSIZE][BRAMSIZE],
	UDTYPE modulus, UDTYPE const_ratio[2])
{
#pragma HLS INLINE off

	UDTYPE operand_temp1, operand_temp2;
	UDTYPE result_temp;

	coff_loop_square:
	for(size_t j = 0; j < BRAMSIZE ;j++){
#pragma HLS PIPELINE
		for (size_t i = 0; i < BRAMNUM; i++){
			operand_temp1 = operand1[i][j];
			operand_temp2 = operand2[i][j];
			result_temp =  HMult_single_return(operand_temp1, operand_temp2, modulus, const_ratio);
			result[i][j] = Hadd_single(result_temp, result_temp, modulus);
		}
	}
}*/

/*void square_N_all(UDTYPE operand1[BRAMNUM][BRAMSIZE],
		UDTYPE operand2[BRAMNUM][BRAMSIZE],
		UDTYPE result1[BRAMSIZE][BRAMSIZE],
		UDTYPE result12[BRAMSIZE][BRAMSIZE],
		UDTYPE result2[BRAMSIZE][BRAMSIZE],
		UDTYPE modulus, UDTYPE modulus_ratio_0, UDTYPE modulus_ratio_1)
{
#pragma HLS INLINE off

	UDTYPE operand_temp1, operand_temp2;
	UDTYPE result_temp1, result_temp2, result_temp12;

	coff_loop_square:
	for(size_t j = 0; j < BRAMSIZE ;j++){
		for (size_t i = 0; i < BRAMNUM; i++){
#pragma HLS PIPELINE II=1
			operand_temp1 = operand1[i][j];
			operand_temp2 = operand2[i][j];

			result_temp1 =  HMult_single_return(operand_temp1, operand_temp1, modulus, modulus_ratio_0, modulus_ratio_1);
			result_temp2 =  HMult_single_return(operand_temp2, operand_temp2, modulus, modulus_ratio_0, modulus_ratio_1);
			result_temp12 =  HMult_single_return(operand_temp1, operand_temp2, modulus, modulus_ratio_0, modulus_ratio_1);

			result1[i][j] = result_temp1;
			result2[i][j] = result_temp2;
			result12[i][j] = Hadd_single(result_temp12, result_temp12, modulus);
		}
	}
}*/
