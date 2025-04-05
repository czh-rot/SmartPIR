#include "top.h"
#include <string>

using namespace std;

// 16384
void cryptonets(UDTYPEin *data_in0) // UDTYPEin *data_in1, UDTYPEin *data_in2, // 56bit unsigned
				// UDTYPEin *data_out0, UDTYPEin *rp_in0, UDTYPEin *key_in0)
{
//#pragma HLS ALLOCATION instances=square_mc3 limit=1 function
// #pragma HLS INTERFACE s_axilite register port = return
//#pragma HLS INTERFACE m_axi = data_in0 depth = 16384 offset = slave bundle = in1

	//--------Modulus DATA init--------------------------------------------------
	UDTYPE modulus[INITMODCOUNT] = {// 66404353,
									// 66420737,
									// 66551809,
									66813953,
									66961409,
									66994177,
									67043329};

	UDTYPE modcopy[INITMODCOUNT] = {// 66404353,
									// 66420737,
									// 66551809,
									66813953,
									66961409,
									66994177,
									67043329};

	UDTYPE modulus_ratio[INITMODCOUNT][3] = {//{66404353, 11391767, 4},
											 //{66420737, 11124097, 4},
											 //{66551809, 8987482, 4},
											 {66813953, 4739403, 4},
											 {66961409, 2364475, 4},
											 {66994177, 1838133, 4},
											 {67043329, 1049584, 4}};

	// --- buffers setting for data & query
	UDTYPE query[MODCOUNT4][2][BRAMNUM][BRAMSIZE]; // 3*2*8*2 = 96
#pragma HLS ARRAY_PARTITION variable = query complete dim = 1
#pragma HLS ARRAY_PARTITION variable = query complete dim = 2
#pragma HLS ARRAY_PARTITION variable = query complete dim = 3
	UDTYPE data[MODCOUNT4][BRAMNUM][BRAMSIZE]; // 3*8*2 = 48
#pragma HLS ARRAY_PARTITION variable = data complete dim = 1
#pragma HLS ARRAY_PARTITION variable = data complete dim = 2
	UDTYPE T[MODCOUNT4][2][BRAMNUM][BRAMSIZE]; // 3*2*8*2 = 96
#pragma HLS ARRAY_PARTITION variable = query complete dim = 1
#pragma HLS ARRAY_PARTITION variable = query complete dim = 2
#pragma HLS ARRAY_PARTITION variable = query complete dim = 3

	// --- input data
	for (IDXTYPE idx = 0; idx < NN; idx++)
	{
#pragma HLS PIPELINE
		IDXTYPE i, j, k;
		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);

		// 0~8191 = 0
		// 8192~16383 = 1
		k = idx >> L_N;

		UDTYPEin in_temp;
		// 0 ~ 16384
		in_temp = data_in0[idx];

		// input for 3 small modulus
		// [3][2][8][1024]
		query[0][k][i][j] = in_temp((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH); // 0-27
		query[1][k][i][j] = in_temp((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH); // 28-55
		query[2][k][i][j] = in_temp((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH); // 56-83
	}

	for (IDXTYPE idx = 0; idx < N; idx++)
	{
#pragma HLS PIPELINE
		IDXTYPE i, j;
		i = idx(L_BRAMNUM - 1, 0);
		j = (idx >> L_BRAMNUM)(L_BRAMSIZE - 1, 0);

		UDTYPEin in_temp;
		// 16384 ~ 24575
		in_temp = data_in0[idx + NN];

		// input for 3 small modulus
		// [3][8][1024]
		data[0][i][j] = in_temp((0 + 1) * BITWIDTH - 1, 0 * BITWIDTH);
		data[1][i][j] = in_temp((1 + 1) * BITWIDTH - 1, 1 * BITWIDTH);
		data[2][i][j] = in_temp((2 + 1) * BITWIDTH - 1, 2 * BITWIDTH);
	}

	// data_test
	// IDXTYPE idx_data = 0;
	pcmult_inplace<MODCOUNT4,MODCOUNT4,BRAMNUM,BRAMSIZE>(query, data, modulus);
	add_inplace(T, query, query);
}
