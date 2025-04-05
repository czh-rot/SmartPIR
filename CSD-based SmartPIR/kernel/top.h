#pragma once

#include <fstream>
#include <iostream>
#include "define.h"
#include "components.h"
#include "multiply.h"
#include "ntt.h"
#include "intt.h"
#include "rescale.h"
#include "add.h"
#include "switchkey.h"
#include "rotation.h"

/*void cryptonets(UDTYPEin *data_in0, UDTYPEin *data_in1, UDTYPEin *data_in2,
		UDTYPEin *data_out0, UDTYPEin *outp,
		UDTYPEin *k0p);*/
void cryptonets(UDTYPEin *data_in0, // UDTYPEin *data_in1, UDTYPEin *data_in2,
				UDTYPEin *data_out0, UDTYPEin *rp_in0, UDTYPEin *key_in0);

template <unsigned modcount>
void init_data_polypoly(UDTYPE operand[modcount][2][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	for (int j = 0; j < BRAMSIZE; j++)
	{
#pragma HLS PIPELINE
		for (int i = 0; i < BRAMNUM; i++)
		{
			for (int m = 0; m < modcount; m++)
			{
				for (int k = 0; k < 2; k++)
				{

					operand[m][k][i][j] = i;
				}
			}
		}
	}
}

template <unsigned modcount>
void init_data_poly(UDTYPE operand[modcount][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	for (int j = 0; j < BRAMSIZE; j++)
	{
#pragma HLS PIPELINE
		for (int i = 0; i < BRAMNUM; i++)
		{
			for (int m = 0; m < modcount; m++)
			{

				operand[m][i][j] = 0;
			}
		}
	}
}

template <unsigned modcount, unsigned polycount>
void data_copy(UDTYPE operand[modcount][polycount][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	for (int m = 0; m < modcount; m++)
	{
		for (int k = 0; k < polycount; k++)
		{
			for (int j = 0; j < BRAMSIZE; j++)
			{
				for (int i = 0; i < BRAMNUM; i++)
				{
#pragma HLS PIPELINE
					operand[m][k][i][j] = m + k + i + j;
				}
			}
		}
	}
}

template <unsigned modcount>
void data_copy_poly1(UDTYPE operand[modcount][BRAMNUM][BRAMSIZE])
{
#pragma HLS INLINE off

	for (int m = 0; m < modcount; m++)
	{
		for (int j = 0; j < BRAMSIZE; j++)
		{
			for (int i = 0; i < BRAMNUM; i++)
			{
#pragma HLS PIPELINE
				operand[m][i][j] = m + i + j;
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------

/*
template <unsigned modcount, unsigned bramnum, unsigned bramsize>
void buffer_copy_polypoly(UDTYPE operand_in[modcount][2][bramnum][bramsize],
		UDTYPE operand_to[modcount][2][bramnum][bramsize]){

	copy:
	for (int j = 0; j < BRAMSIZE; j++){
#pragma HLS PIPELINE
		for(int i = 0; i < BRAMNUM; i++){
			for (int k = 0; k < 2; k++){
				for(int m = 0; m < modcount; m++){
					b2[3][k][i][j] = b1_o[4][k][i][j];
					b2[2][k][i][j] = b1_o[3][k][i][j];
					b2[1][k][i][j] = b1_o[2][k][i][j];
					b1[2][k][i][j] = b1_o[1][k][i][j];
					b1[1][k][i][j] = b1_o[0][k][i][j];

				}

			}
		}
	}
}
*/

// 16384 cycles
template <unsigned bramnum, unsigned bramsize>
void layer3_out_mc2(UDTYPE operand0[2][bramnum][bramsize],
					UDTYPE operand1[2][bramnum][bramsize],
					UDTYPEin *out)
{
#pragma HLS INLINE off

	UDTYPEin temp_send;

	IDXTYPE i, j, k;
	for (IDXTYPE idx = 0; idx < NN; idx++)
	{
#pragma HLS PIPELINE

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);
		k = idx >> L_N;

		temp_send((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH) = operand0[k][i][j];
		temp_send((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH) = operand1[k][i][j];

		out[idx] = temp_send;
	}
}

// 16384 cycles
template <unsigned bramnum, unsigned bramsize>
void layer5_out_mc1(UDTYPE operand0[2][bramnum][bramsize], UDTYPEin *out)
{
#pragma HLS INLINE off

	UDTYPEin temp_send;

	/*for (int k = 0; k < 2; k++){
		for (int j = 0; j < bramsize; j++){
			for (int i = 0; i < bramnum; i++){
#pragma HLS PIPELINE

				temp_send((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH) = operand0[k][i][j];

				*out = temp_send;
				out++;
			}
		}
	}*/
	IDXTYPE i, j, k;
	for (IDXTYPE idx = 0; idx < NN; idx++)
	{
#pragma HLS PIPELINE

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);
		k = idx >> L_N;

		temp_send((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH) = operand0[k][i][j];

		out[idx] = temp_send;
	}
}

template <unsigned bramnum, unsigned bramsize>
void layer4_datain_mc3(UDTYPEin *in_k0, UDTYPEin *in_k1, UDTYPE operand[3][2][bramnum][bramsize])
{
#pragma HLS INLINE off
	UDTYPEin in_temp0, in_temp1;
	UDTYPEin op_tempk0_m0, op_tempk0_m1, op_tempk0_m2;
	UDTYPEin op_tempk1_m0, op_tempk1_m1, op_tempk1_m2;

	/*int index = 0;
	copy_data_in:
	for (size_t j = 0; j < bramsize; j++){
		for (size_t i = 0; i < bramnum; i++){
#pragma HLS PIPELINE
			in_temp0 = in_k0[index];
			in_temp1 = in_k1[index];

			op_tempk0_m0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
			op_tempk0_m1 = (in_temp0)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
			op_tempk0_m2 = (in_temp0)((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);

			operand[0][0][i][j] = op_tempk0_m0;
			operand[1][0][i][j] = op_tempk0_m1;
			operand[2][0][i][j] = op_tempk0_m2;

			op_tempk1_m0 = (in_temp1)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
			op_tempk1_m1 = (in_temp1)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
			op_tempk1_m2 = (in_temp1)((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);

			operand[0][1][i][j] = op_tempk1_m0;
			operand[1][1][i][j] = op_tempk1_m1;
			operand[2][1][i][j] = op_tempk1_m2;

			index++;
		}
	}*/

	IDXTYPE i, j;
copy_data_in:
	for (IDXTYPE idx = 0; idx < N; idx++)
	{
#pragma HLS PIPELINE

		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);

		in_temp0 = in_k0[idx];
		in_temp1 = in_k1[idx];

		op_tempk0_m0 = (in_temp0)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		op_tempk0_m1 = (in_temp0)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		op_tempk0_m2 = (in_temp0)((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);

		operand[0][0][i][j] = op_tempk0_m0;
		operand[1][0][i][j] = op_tempk0_m1;
		operand[2][0][i][j] = op_tempk0_m2;

		op_tempk1_m0 = (in_temp1)((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		op_tempk1_m1 = (in_temp1)((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		op_tempk1_m2 = (in_temp1)((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);

		operand[0][1][i][j] = op_tempk1_m0;
		operand[1][1][i][j] = op_tempk1_m1;
		operand[2][1][i][j] = op_tempk1_m2;
	}
}

template <unsigned modcount>
void root_powers_copy(UDTYPE rp[modcount][N], UDTYPE srp[modcount][N],
					  UDTYPE idt[modcount][N], UDTYPE sidt[modcount][N])
{
#pragma HLS INLINE off

root_copy_i:
	for (int i = 0; i < modcount; i++)
	{
	root_copy_j:
		for (int j = 0; j < N; j++)
		{
#pragma HLS PIPELINE
			rp[i][j] = i + j;
			srp[i][j] = i + j;
			idt[i][j] = i + j;
			sidt[i][j] = i + j;
		}
	}
}

template <unsigned bramnum, unsigned bramsize>
void send_data_N(UDTYPE from[BRAMNUM][BRAMSIZE], UDTYPEin *to)
{
#pragma HLS INLINE off
	UDTYPEin temp_send;
	for (uint32_t j = 0; j < BRAMSIZE; j++)
	{
		for (uint32_t i = 0; i < BRAMNUM; i += 2)
		{
#pragma HLS PIPELINE

			temp_send((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH) = from[i][j / BRAMNUM];
			i++;
			temp_send((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH) = from[i][j / BRAMNUM];
			*to = temp_send;
			to++;
		}
	}
}

template <unsigned modcount, unsigned bramnum, unsigned bramsize>
void send_data_polys(UDTYPE from[modcount][2][bramnum][bramsize], UDTYPEin *to)
{
#pragma HLS INLINE off
axi_send:
	for (uint32_t m = 0; m < modcount; m++)
	{
		for (uint32_t k = 0; k < 2; k++)
		{
			send_data_N<bramnum, bramsize>(from[m][k], to);
		}
	}
}

/*template<unsigned from_size>
void copy_data_1d(UDTYPE *from,UDTYPE *to){

	for(int i=0;i<from_size;i++){
		to[i]=from[i];
	}
}

template<unsigned from_size>
void copy_data_1d_1to2(UDTYPEin *from,UDTYPE *to){
	for(int i=0;i<from_size;i = i + 2){
		UDTYPE2 tmp;
		tmp = *from++;
		to[i]=tmp(BITWIDTH-1,0);
		to[i + 1]=tmp(2*BITWIDTH-1,BITWIDTH);
	}
}

template<unsigned i_size,unsigned j_size>
void copy_data_2d_1to2(UDTYPEin *from,UDTYPE to1[i_size][j_size],UDTYPE to2[i_size][j_size]){
	for(int i=0;i<i_size;i++){
		for(int j=0;j<j_size;j = j + 2){
#pragma HLS PIPELINE
			to1[i][j]=(*from)(BITWIDTH-1,0);
			to1[i][j + 1]=(*from)(2*BITWIDTH-1,BITWIDTH);
			to2[i][j]=(*from)(BITWIDTH-1,0);
			to2[i][j + 1]=(*from)(2*BITWIDTH-1,BITWIDTH);
			from++;
		}
	}
}
template<unsigned i_size,unsigned j_size>
void copy_data_2d_1to2(UDTYPEin *from,UDTYPE to1[i_size][j_size]){
	for(int i=0;i<i_size;i++){
		for(int j=0;j<j_size;j = j + 2){
#pragma HLS PIPELINE
			to1[i][j]=(*from)(BITWIDTH-1,0);
			to1[i][j + 1]=(*from)(2*BITWIDTH-1,BITWIDTH);
			from++;
		}
	}
}

template<unsigned i_size,unsigned j_size>
void copy_data_2d(UDTYPEin *from,UDTYPE to1[i_size][j_size]){
	for(int i=0;i<i_size;i++){
		for(int j=0;j<j_size;j = j ++){
#pragma HLS PIPELINE
			to1[i][j]=(*from)(BITWIDTH-1,0);
			//to1[i][j + 1]=(*from)(2*BITWIDTH-1,BITWIDTH);
			from++;
		}
	}
}*/

// template<unsigned i_size,unsigned j_size,unsigned k_size>
// void copy_data_3d(UDTYPE *from,UDTYPE to[i_size][j_size][k_size]){
//	for(int i=0;i<i_size;i++){
//		for(int j=0;j<j_size;j++){
//			for(int k=0;k<k_size;k++){
//				to[i][j][k]=*from++;
//			}
//		}
//	}
// }
/*template<unsigned i_size,unsigned j_size,unsigned k_size>
void copy_data_3d(UDTYPEin *from,UDTYPE to[i_size][j_size][k_size]){
	for(int i=0;i<i_size;i++){
		for(int j=0;j<j_size;j++){
			for(int k=0;k<k_size;k++){
				to[i][j][k]=(*from)(BITWIDTH-1,0);
				from++;
			}
		}
	}
}

template<unsigned i_size,unsigned j_size,unsigned k_size,unsigned l_size>
void copy_data_4d(UDTYPEin *from,UDTYPE to[i_size][j_size][k_size][l_size]){
	for(int i=0;i<i_size;i++){
		for(int j=0;j<j_size;j++){
			for(int k=0;k<k_size;k++){
				for(int l=0;l<l_size;l++){
				to[i][j][k][l]=(*from)(BITWIDTH-1,0);
				from++;
				}
			}
		}
	}
}
template<unsigned i_size,unsigned j_size>
void copy_data_2d(UDTYPE from[i_size][j_size],UDTYPE to[i_size][j_size]){

for(int j=0;j<j_size;j++){
#pragma HLS PIPELINE
			for(int i=0;i<i_size;i++){

				to[i][j]=from[i][j];

		}
	}
}

template<unsigned i_size,unsigned j_size,unsigned k_size>
void copy_data_3d_1to2(UDTYPEin *from,UDTYPE to[i_size][j_size][k_size]){
	for(int i=0;i<i_size;i++){
		for(int j=0;j<j_size;j++){
			for(int k=0;k<k_size;k = k + 2){
				to[i][j][k]=(*from)(BITWIDTH-1,0);
				to[i][j][k + 1]=(*from)(2*BITWIDTH-1,BITWIDTH);
				from++;
			}
		}
	}
}

template<unsigned i_size,unsigned j_size,unsigned k_size,unsigned l_size>
void copy_data_4d_1to2(UDTYPEin *from,UDTYPE to[i_size][j_size][k_size][l_size]){
	for(int i=0;i<i_size;i++){
		for(int j=0;j<j_size;j++){
			for(int k=0;k<k_size;k++){
				for(int l=0;l<l_size;l+=2){
					to[i][j][k][l]=(*from)(BITWIDTH-1,0);
					to[i][j][k][l+1]=(*from)(2*BITWIDTH-1,BITWIDTH);
					from++;
				}
			}
		}
	}
}




template <unsigned mod_count, unsigned relinsize>
void copy_data_block_conv(UDTYPEin *from3,UDTYPE to1[mod_count][N],UDTYPE to2[mod_count][N],UDTYPE to3[mod_count][N],UDTYPE to4[mod_count][N],
		UDTYPE to5[mod_count][N],UDTYPE to6[mod_count][N],UDTYPE to7[mod_count][N],UDTYPE to8[mod_count][N], UDTYPE to9[relinsize]){

	//copy_data_3d_1to2<pic_size,polys_num,poly_size>(from1,to1);
	//copy_data_2d_1to2<pic_size,poly_size>(from2,to2);

	copy_data_2d_1to2<mod_count,N>(from3,to1,to2);
	copy_data_2d_1to2<mod_count,N>(from3 + mod_count * (N / 2), to3, to4);
	copy_data_2d_1to2<mod_count,N>(from3 + mod_count * (N / 2 * 2), to5, to6);
	copy_data_2d_1to2<mod_count,N>(from3 + mod_count * (N / 2 * 3), to7, to8);
	copy_data_1d_1to2<RELINKEYS_SIZE>(from3 + mod_count * (N / 2 * 4), to9);

}

template <unsigned mod_count>
void copy_data_block_conv(UDTYPEin *from3,UDTYPE to1[mod_count][N],UDTYPE to2[mod_count][N],UDTYPE to3[mod_count][N],UDTYPE to4[mod_count][N],
		UDTYPE to5[mod_count][N],UDTYPE to6[mod_count][N],UDTYPE to7[mod_count][N],UDTYPE to8[mod_count][N]){

	//copy_data_3d_1to2<pic_size,polys_num,poly_size>(from1,to1);
	//copy_data_2d_1to2<pic_size,poly_size>(from2,to2);

	copy_data_2d_1to2<mod_count,N>(from3,to1,to2);
	copy_data_2d_1to2<mod_count,N>(from3 + mod_count * (N / 2), to3, to4);
	copy_data_2d_1to2<mod_count,N>(from3 + mod_count * (N / 2 * 2), to5, to6);
	copy_data_2d_1to2<mod_count,N>(from3 + mod_count * (N / 2 * 3), to7, to8);
	//copy_data_1d_1to2<RELINKEYS_SIZE>(from3 + mod_count * (N / 2 * 4), to9);

}

template<unsigned mod_count>
void copy_data_switchkey(UDTYPEin *from1,UDTYPEin * from2,UDTYPE to1[mod_count][N],UDTYPE to2[mod_count][N],UDTYPE to3[mod_count][N],UDTYPE to4[mod_count][N],UDTYPE to5[MODCOUNT2][2][BRAMNUM][BRAMSIZE],UDTYPE to6[MODCOUNT2][BRAMNUM][BRAMSIZE]){

	//copy_data_3d_1to2<pic_size,polys_num,poly_size>(from1,to1);
	//copy_data_2d_1to2<pic_size,poly_size>(from2,to2);

	copy_data_2d_1to2<mod_count,N>(from1,to1);
	copy_data_2d_1to2<mod_count,N>(from1+ mod_count * (N / 2),to2);
	copy_data_2d_1to2<mod_count,N>(from1+ mod_count * (N / 2 * 2),to3);
	copy_data_2d_1to2<mod_count,N>(from1+ mod_count * (N / 2 * 3),to4);

	int index=81920;

	copy_data_4d<mod_count,2,BRAMNUM,BRAMSIZE>(from2,to5);

	copy_data_3d<mod_count,BRAMNUM,BRAMSIZE>(from2+index,to6);

}
template<unsigned init_mod_count,unsigned mod_count>
void copy_data_switchkey_layer3(UDTYPEin *from1,UDTYPEin * from2,UDTYPE to1[init_mod_count][N],UDTYPE to2[init_mod_count][N],UDTYPE to3[init_mod_count][N],UDTYPE to4[init_mod_count][N],UDTYPE to5[mod_count][2][BRAMNUM][BRAMSIZE],UDTYPE to6[mod_count][BRAMNUM][BRAMSIZE]){

	//copy_data_3d_1to2<pic_size,polys_num,poly_size>(from1,to1);
	//copy_data_2d_1to2<pic_size,poly_size>(from2,to2);

	copy_data_2d_1to2<init_mod_count,N>(from1,to1);
	copy_data_2d_1to2<init_mod_count,N>(from1+ mod_count * (N / 2),to2);
	copy_data_2d_1to2<init_mod_count,N>(from1+ mod_count * (N / 2 * 2),to3);
	copy_data_2d_1to2<init_mod_count,N>(from1+ mod_count * (N / 2 * 3),to4);

	//int index=81920;
	int index=8192*8;

	copy_data_4d<mod_count,2,BRAMNUM,BRAMSIZE>(from2,to5);

	copy_data_3d<mod_count,BRAMNUM,BRAMSIZE>(from2+index,to6);

}

template <unsigned polys_num,unsigned mod_count>  //2,6
void copy_data_encrypted(UDTYPEin *from,UDTYPE to[polys_num][mod_count][BRAMNUM][BRAMSIZE]){

	for(size_t polys_index = 0; polys_index < polys_num; polys_index++){
		for(size_t mod_index = 0; mod_index < mod_count; mod_index++){
			copy_data_N(from,to[polys_index][mod_index][0],to[polys_index][mod_index][1],to[polys_index][mod_index][2],to[polys_index][mod_index][3],
					to[polys_index][mod_index][4],to[polys_index][mod_index][5],to[polys_index][mod_index][6],to[polys_index][mod_index][7]);
			from = from + N/2;
		}
	}
}


template <unsigned polys_num,unsigned mod_count>
void send_data_encrypted(UDTYPE from[polys_num][mod_count][BRAMNUM][BRAMSIZE], UDTYPEin *to){

	for(size_t polys_index = 0; polys_index < polys_num; polys_index++){
		for(size_t mod_index = 0; mod_index < mod_count; mod_index++){
			send_data_N(from[polys_index][mod_index],to);
			to = to + N/2;
		}
	}
}*/
