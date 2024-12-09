#include "define.h"
using namespace std;

void inverse_ntt_negacyclic_harvey_lazy_core(UDTYPE W, UDTYPE Wprime, UDTYPE modulus, UDTYPE two_times_modulus,
											 UDTYPE &operanda, UDTYPE &operandb);

void inverse_ntt_core(UDTYPE W, UDTYPE Wprime, UDTYPE modulus, UDTYPE two_times_modulus,
					  UDTYPE operanda_in, UDTYPE operandb_in, UDTYPE &operanda, UDTYPE &operandb);

template <unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax>
void intt_1core(UDTYPE operand[bramnum][bramsize],
				UDTYPE idt[N], UDTYPE sidt[N], UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				uint32_t x = (arrayWindex[m] + temp_adder);
				//				w[m]=rp[x];
				//				sw[m]=srp[x];
				w[m] = idt[0];
				sw[m] = sidt[0];
			}

			operanda[0] = operand[0][MEa];
			operandb[0] = operand[1][MEb];

			inverse_ntt_negacyclic_harvey_lazy_core(w[0], sw[0], modulus, two_times_modulus, operanda[0], operandb[0]);

			operand[0][MEa] = operanda[0];
			operand[1][MEb] = operandb[0];

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1:
	for (uint32_t i = stage_num - stage1_max; i < stage_num; i++)
	{
		uint32_t stepsize = step_num >> (i + 1);

		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			UDTYPE w0 = idt[temp1 + (j >> temp2)];
			UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			operanda[0] = operand[0][MEa];
			operandb[0] = operand[0][MEb];
			operanda_[0] = operand[0 + corenum][MEa];
			operandb_[0] = operand[0 + corenum][MEb];

			inverse_ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[0], operandb[0]);

			j++;

			inverse_ntt_negacyclic_harvey_lazy_core(w00, sw00, modulus, two_times_modulus, operanda_[0], operandb_[0]);

			operand[0][MEa] = operanda[0];
			operand[0][MEb] = operandb[0];
			operand[0 + corenum][MEa] = operanda_[0];
			operand[0 + corenum][MEb] = operandb_[0];

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}
}

template <unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax>
void intt_2core(UDTYPE operand[bramnum][bramsize],
				UDTYPE idt[N], UDTYPE sidt[N], UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				uint32_t x = (arrayWindex[m] + temp_adder);
				//				w[m]=rp[x];
				//				sw[m]=srp[x];
				w[m] = idt[0];
				sw[m] = sidt[0];
			}

			if (i == 0)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[1][MEb];
				operanda[1] = operand[2][MEa];
				operandb[1] = operand[3][MEb];
			}
			else if (i == 1)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[2][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[3][MEb];
			}

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w[m], sw[m], modulus, two_times_modulus, operanda[m], operandb[m]);
			}

			if (i == 0)
			{
				operand[0][MEa] = operanda[0];
				operand[1][MEb] = operandb[0];
				operand[2][MEa] = operanda[1];
				operand[3][MEb] = operandb[1];
			}
			else if (i == 1)
			{
				operand[0][MEa] = operanda[0];
				operand[2][MEb] = operandb[0];
				operand[1][MEa] = operanda[1];
				operand[3][MEb] = operandb[1];
			}

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1:
	for (uint32_t i = stage_num - stage1_max; i < stage_num; i++)
	{
		uint32_t stepsize = step_num >> (i + 1);

		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			UDTYPE w0 = idt[temp1 + (j >> temp2)];
			UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operanda[m] = operand[m][MEa];
				operandb[m] = operand[m][MEb];
				operanda_[m] = operand[m + corenum][MEa];
				operandb_[m] = operand[m + corenum][MEb];
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m]);
			}
			j++;
			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m]);
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operand[m][MEa] = operanda[m];
				operand[m][MEb] = operandb[m];
				operand[m + corenum][MEa] = operanda_[m];
				operand[m + corenum][MEb] = operandb_[m];
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}
}

template <unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax>
void intt_4core(UDTYPE operand[bramnum][bramsize],
				UDTYPE idt[N], UDTYPE sidt[N], UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE operanda_out[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_out complete dim = 0
	UDTYPE operandb_out[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_out complete dim = 0
	UDTYPE operanda_2out[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_2out complete dim = 0
	UDTYPE operandb_2out[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_2out complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				uint32_t x = (arrayWindex[m] + temp_adder);
				//				w[m]=rp[x];
				//				sw[m]=srp[x];
				w[m] = idt[0];
				sw[m] = sidt[0];
			}

			if (i == 0)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[1][MEb];
				operanda[1] = operand[2][MEa];
				operandb[1] = operand[3][MEb];
				operanda[2] = operand[4][MEa];
				operandb[2] = operand[5][MEb];
				operanda[3] = operand[6][MEa];
				operandb[3] = operand[7][MEb];

				for (uint32_t m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus, operanda[m], operandb[m], operanda_[m], operandb_[m]);
				}

				operand[0][MEa] = operanda_[0];
				operand[1][MEb] = operandb_[0];
				operand[2][MEa] = operanda_[1];
				operand[3][MEb] = operandb_[1];
				operand[4][MEa] = operanda_[2];
				operand[5][MEb] = operandb_[2];
				operand[6][MEa] = operanda_[3];
				operand[7][MEb] = operandb_[3];
			}
			else if (i == 1)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[2][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[3][MEb];
				operanda[2] = operand[4][MEa];
				operandb[2] = operand[6][MEb];
				operanda[3] = operand[5][MEa];
				operandb[3] = operand[7][MEb];

				for (uint32_t m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus, operanda[m], operandb[m], operanda_[m], operandb_[m]);
				}

				operand[0][MEa] = operanda_[0];
				operand[2][MEb] = operandb_[0];
				operand[1][MEa] = operanda_[1];
				operand[3][MEb] = operandb_[1];
				operand[4][MEa] = operanda_[2];
				operand[6][MEb] = operandb_[2];
				operand[5][MEa] = operanda_[3];
				operand[7][MEb] = operandb_[3];
			}
			else if (i == 2)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[4][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[5][MEb];
				operanda[2] = operand[2][MEa];
				operandb[2] = operand[6][MEb];
				operanda[3] = operand[3][MEa];
				operandb[3] = operand[7][MEb];

				for (uint32_t m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus, operanda[m], operandb[m], operanda_[m], operandb_[m]);
				}

				operand[0][MEa] = operanda_[0];
				operand[4][MEb] = operandb_[0];
				operand[1][MEa] = operanda_[1];
				operand[5][MEb] = operandb_[1];
				operand[2][MEa] = operanda_[2];
				operand[6][MEb] = operandb_[2];
				operand[3][MEa] = operanda_[3];
				operand[7][MEb] = operandb_[3];
			}

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1:
	for (uint32_t i = stage_num - stage1_max; i < stage_num; i++)
	{
		// uint32_t stepsize = step_num >> (i + 1);
		uint32_t stepsize = 1 << (i - (stage_num - stage1_max));

		/*std::cout << " stepsize = " <<stepsize << std::endl;
		std::cout << " i = " <<i << std::endl;
		std::cout << " stage_num = " <<stage_num << std::endl;
		std::cout << " stage1_max = " <<stage1_max << std::endl;*/

		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			UDTYPE w0 = idt[temp1 + (j >> temp2)];
			UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operanda[m] = operand[m][MEa];
				operandb[m] = operand[m][MEb];
				operanda_[m] = operand[m + corenum][MEa];
				operandb_[m] = operand[m + corenum][MEb];
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL

				inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m],
								 operanda_out[m], operandb_out[m]);

				inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m],
								 operanda_2out[m], operandb_2out[m]);
				// inverse_ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m]);
			}
			j++;
			/*for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m]);
			}*/

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL

				/*std::cout << " m = " <<m << std::endl;
				std::cout << " m + corenum = " <<m + corenum << std::endl;
				std::cout << " MEa = " <<MEa << std::endl;
				std::cout << " MEb = " <<MEb << std::endl;*/

				operand[m][MEa] = operanda_out[m];
				operand[m][MEb] = operandb_out[m];
				operand[m + corenum][MEa] = operanda_2out[m];
				operand[m + corenum][MEb] = operandb_2out[m];
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}
}

template <unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax, unsigned type>
void intt_4core_new(UDTYPE operand[bramnum][bramsize],
					UDTYPE idt[RPBRAMNUM][RPBRAMSIZE],
					UDTYPE sidt[RPBRAMNUM][RPBRAMSIZE],
					UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE operanda_out[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_out complete dim = 0
	UDTYPE operandb_out[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_out complete dim = 0
	UDTYPE operanda_2out[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_2out complete dim = 0
	UDTYPE operandb_2out[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_2out complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			// 			for (uint32_t m = 0; m < corenum; m++)
			// 			{
			// #pragma HLS UNROLL
			// 				uint32_t x = (arrayWindex[m] + temp_adder);
			// 				//				w[m]=rp[x];
			// 				//				sw[m]=srp[x];
			// 				w[m] = idt[0];
			// 				sw[m] = sidt[0];
			// 			}

			if (i == 0)
			{
				IDXTYPE x = (arrayWindex[0] + temp_adder);
				IDXTYPE rp_idx0_j;
				rp_idx0_j = x >> L_RPBRAMNUM;

				w[0] = idt[0][rp_idx0_j];
				sw[0] = sidt[0][rp_idx0_j];
				w[1] = idt[1][rp_idx0_j];
				sw[1] = sidt[1][rp_idx0_j];
				w[2] = idt[2][rp_idx0_j];
				sw[2] = sidt[2][rp_idx0_j];
				w[3] = idt[3][rp_idx0_j];
				sw[3] = sidt[3][rp_idx0_j];

				operanda[0] = operand[0][MEa];
				operandb[0] = operand[1][MEb];
				operanda[1] = operand[2][MEa];
				operandb[1] = operand[3][MEb];
				operanda[2] = operand[4][MEa];
				operandb[2] = operand[5][MEb];
				operanda[3] = operand[6][MEa];
				operandb[3] = operand[7][MEb];

				for (uint32_t m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL

					inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
									 operanda[m], operandb[m], operanda_[m], operandb_[m]);
				}

				operand[0][MEa] = operanda_[0];
				operand[1][MEb] = operandb_[0];
				operand[2][MEa] = operanda_[1];
				operand[3][MEb] = operandb_[1];
				operand[4][MEa] = operanda_[2];
				operand[5][MEb] = operandb_[2];
				operand[6][MEa] = operanda_[3];
				operand[7][MEb] = operandb_[3];
			}
			else if (i == 1)
			{
				IDXTYPE x = (arrayWindex[0] + temp_adder);
				IDXTYPE rp_idx0_i, rp_idx0_j;
				rp_idx0_i = x(L_RPBRAMNUM - 1, 0);
				rp_idx0_j = x >> L_RPBRAMNUM;

				UDTYPE rp_t1, srp_t1;
				UDTYPE rp_t2, srp_t2;

				rp_t1 = idt[rp_idx0_i][rp_idx0_j];
				srp_t1 = sidt[rp_idx0_i][rp_idx0_j];

				rp_t2 = idt[rp_idx0_i + 1][rp_idx0_j];
				srp_t2 = sidt[rp_idx0_i + 1][rp_idx0_j];

				w[0] = rp_t1;
				sw[0] = srp_t1;
				w[1] = rp_t1;
				sw[1] = srp_t1;
				w[2] = rp_t2;
				sw[2] = srp_t2;
				w[3] = rp_t2;
				sw[3] = srp_t2;

				operanda[0] = operand[0][MEa];
				operandb[0] = operand[2][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[3][MEb];
				operanda[2] = operand[4][MEa];
				operandb[2] = operand[6][MEb];
				operanda[3] = operand[5][MEa];
				operandb[3] = operand[7][MEb];

				for (uint32_t m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					//					cout << "w[m]" << w[m] << endl;
					//					cout << "sw[m]" << sw[m] << endl;
					//
					//					cout << "in: operanda[m]" << operanda[m] << endl;
					//					cout << "in: operandb[m]" << operandb[m] << endl;
					inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
									 operanda[m], operandb[m], operanda_[m], operandb_[m]);

					//					cout << "out: operanda[m]" << operanda_[m] << endl;
					//					cout << "out: operandb[m]" << operandb_[m] << endl;
				}

				operand[0][MEa] = operanda_[0];
				operand[2][MEb] = operandb_[0];
				operand[1][MEa] = operanda_[1];
				operand[3][MEb] = operandb_[1];
				operand[4][MEa] = operanda_[2];
				operand[6][MEb] = operandb_[2];
				operand[5][MEa] = operanda_[3];
				operand[7][MEb] = operandb_[3];
			}
			else if (i == 2)
			{
				for (uint32_t m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL

					IDXTYPE x = (arrayWindex[0] + temp_adder);
					IDXTYPE rp_idx0_i, rp_idx0_j;
					rp_idx0_i = x(L_RPBRAMNUM - 1, 0);
					rp_idx0_j = x >> L_RPBRAMNUM;
					// cout << "x = " << x << endl;

					w[m] = idt[rp_idx0_i][rp_idx0_j];
					sw[m] = sidt[rp_idx0_i][rp_idx0_j];
				}

				operanda[0] = operand[0][MEa];
				operandb[0] = operand[4][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[5][MEb];
				operanda[2] = operand[2][MEa];
				operandb[2] = operand[6][MEb];
				operanda[3] = operand[3][MEa];
				operandb[3] = operand[7][MEb];

				for (uint32_t m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					//					cout << "w[m]" << w[m] << endl;
					//					cout << "sw[m]" << sw[m] << endl;
					//
					//					cout << "in: operanda[m]" << operanda[m] << endl;
					//					cout << "in: operandb[m]" << operandb[m] << endl;
					inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
									 operanda[m], operandb[m], operanda_[m], operandb_[m]);
					//					cout << "out: operanda[m]" << operanda_[m] << endl;
					//					cout << "out: operandb[m]" << operandb_[m] << endl;
				}

				operand[0][MEa] = operanda_[0];
				operand[4][MEb] = operandb_[0];
				operand[1][MEa] = operanda_[1];
				operand[5][MEb] = operandb_[1];
				operand[2][MEa] = operanda_[2];
				operand[6][MEb] = operandb_[2];
				operand[3][MEa] = operanda_[3];
				operand[7][MEb] = operandb_[3];
			}

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1_without_last:
	for (uint32_t i = stage_num - stage1_max; i < stage_num - 1; i++)
	{
		// uint32_t stepsize = step_num >> (i + 1);

		// stepsize formula is wrong
		uint32_t stepsize = 1 << (i - (stage_num - stage1_max));

		/*std::cout << " stepsize = " <<stepsize << std::endl;
		std::cout << " i = " <<i << std::endl;
		std::cout << " stage_num = " <<stage_num << std::endl;
		std::cout << " stage1_max = " <<stage1_max << std::endl;*/

		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;
		uint32_t stage_tempa = 12;

		MEa = 0;
		MEb = stepsize;

		// cout << "stepsize" << stepsize << endl;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			// UDTYPE w0 = idt[temp1 + (j >> temp2)];
			// UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			// UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			// UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			IDXTYPE rp_idx0 = temp1 + (j >> temp2);
			IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
			IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

			rp_idx0_i = rp_idx0(L_RPBRAMNUM - 1, 0);
			rp_idx0_j = rp_idx0 >> L_RPBRAMNUM;
			rp_idx00_i = rp_idx00(L_RPBRAMNUM - 1, 0);
			rp_idx00_j = rp_idx00 >> L_RPBRAMNUM;

			UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
			UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
			UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
			UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				//				cout << "m = " << m << endl;
				//				cout << "MEa" << MEa << endl;
				//				cout << "MEb = " << MEb << endl;
				//				cout << "[m][MEa] = " << operand[m][MEa] << endl;
				//				cout << "[m][MEb] = " << operand[m][MEb] << endl;

				operanda[m] = operand[m][MEa];
				operandb[m] = operand[m][MEb];
				operanda_[m] = operand[m + corenum][MEa];
				operandb_[m] = operand[m + corenum][MEb];
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				/*
								if (i == stage_tempa){
									cout << "w[m]" << w0 << endl;
									cout << "sw[m]" << sw0 << endl;

									cout << "in: operanda[m]" << operanda[m] << endl;
									cout << "in: operandb[m]" << operandb[m] << endl;
								}*/

				inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m],
								 operanda_out[m], operandb_out[m]);

				/*if (i == stage_tempa){
				cout << "out: operanda[m]" << operanda_out[m] << endl;
				cout << "out: operandb[m]" << operandb_out[m] << endl;
				}*/

				/*if (i == stage_tempa){
					cout << "w[m]" << w0 << endl;
					cout << "sw[m]" << sw0 << endl;

					cout << "in: operanda[m]" << operanda_[m] << endl;
					cout << "in: operandb[m]" << operandb_[m] << endl;
				}*/

				inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m],
								 operanda_2out[m], operandb_2out[m]);

				/*if (i == stage_tempa){
				cout << "out: operanda[m]" << operanda_2out[m] << endl;
				cout << "out: operandb[m]" << operandb_2out[m] << endl;
				}*/
				// inverse_ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m]);
			}
			j++;

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL

				/*std::cout << " m = " <<m << std::endl;
				std::cout << " m + corenum = " <<m + corenum << std::endl;
				std::cout << " MEa = " <<MEa << std::endl;
				std::cout << " MEb = " <<MEb << std::endl;*/

				operand[m][MEa] = operanda_out[m];
				operand[m][MEb] = operandb_out[m];
				operand[m + corenum][MEa] = operanda_2out[m];
				operand[m + corenum][MEb] = operandb_2out[m];
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}

	if (type == 0)
	{
	type0_intt_stage1_last:
		for (uint32_t i = stage_num - 1; i < stage_num; i++)
		{
			// uint32_t stepsize = step_num >> (i + 1);

			// stepsize formula is wrong
			uint32_t stepsize = 1 << (i - (stage_num - stage1_max));

			/*std::cout << " stepsize = " <<stepsize << std::endl;
			std::cout << " i = " <<i << std::endl;
			std::cout << " stage_num = " <<stage_num << std::endl;
			std::cout << " stage1_max = " <<stage1_max << std::endl;*/

			uint32_t temp2 = i - (stage_num - stage1_max - 1);
			uint32_t temp1 = step_num >> temp2;
			uint32_t temp3 = (one_u64 << temp2) - one_u64;
			uint32_t stage_tempa = 12;

			MEa = 0;
			MEb = stepsize;

			// cout << "stepsize" << stepsize << endl;

		type0_intt_stage1_inner_last:
			for (uint32_t j = 0; j < step_num; j++)
			{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

				// UDTYPE w0 = idt[temp1 + (j >> temp2)];
				// UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
				// UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
				// UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

				IDXTYPE rp_idx0 = temp1 + (j >> temp2);
				IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
				IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

				rp_idx0_i = rp_idx0(L_RPBRAMNUM - 1, 0);
				rp_idx0_j = rp_idx0 >> L_RPBRAMNUM;
				rp_idx00_i = rp_idx00(L_RPBRAMNUM - 1, 0);
				rp_idx00_j = rp_idx00 >> L_RPBRAMNUM;

				UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
				UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
				UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
				UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

				for (int m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					operanda[m] = operand[m][MEa];
					operandb[m] = operand[m][MEb];
					operanda_[m] = operand[m + corenum][MEa];
					operandb_[m] = operand[m + corenum][MEb];
				}

				for (int m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m],
									 operanda_out[m], operandb_out[m]);

					inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m],
									 operanda_2out[m], operandb_2out[m]);
				}
				j++;

				for (int m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					operand[m][MEa] = operanda_out[m];
					operand[m][MEb] = operandb_out[m];
					operand[m + corenum][MEa] = operanda_2out[m];
					operand[m + corenum][MEb] = operandb_2out[m];
				}

				if (((j + 1) & temp3) == 0)
				{
					MEa += stepsize + 1;
					MEb += stepsize + 1;
				}
				else
				{
					MEa += 1;
					MEb += 1;
				}
			}
		}
	}
	else if (type == 1)
	{
	type1_intt_stage1_last:
		for (uint32_t i = stage_num - 1; i < stage_num; i++)
		{
			// uint32_t stepsize = step_num >> (i + 1);

			// stepsize formula is wrong
			uint32_t stepsize = 1 << (i - (stage_num - stage1_max));

			/*std::cout << " stepsize = " <<stepsize << std::endl;
			std::cout << " i = " <<i << std::endl;
			std::cout << " stage_num = " <<stage_num << std::endl;
			std::cout << " stage1_max = " <<stage1_max << std::endl;*/

			uint32_t temp2 = i - (stage_num - stage1_max - 1);
			uint32_t temp1 = step_num >> temp2;
			uint32_t temp3 = (one_u64 << temp2) - one_u64;
			uint32_t stage_tempa = 12;

			MEa = 0;
			MEb = stepsize;

			// cout << "stepsize" << stepsize << endl;

		type1_intt_stage1_inner_last:
			for (uint32_t j = 0; j < step_num; j++)
			{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

				// UDTYPE w0 = idt[temp1 + (j >> temp2)];
				// UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
				// UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
				// UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

				IDXTYPE rp_idx0 = temp1 + (j >> temp2);
				IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
				IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

				rp_idx0_i = rp_idx0(L_RPBRAMNUM - 1, 0);
				rp_idx0_j = rp_idx0 >> L_RPBRAMNUM;
				rp_idx00_i = rp_idx00(L_RPBRAMNUM - 1, 0);
				rp_idx00_j = rp_idx00 >> L_RPBRAMNUM;

				UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
				UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
				UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
				UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

				for (int m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					operanda[m] = operand[m][MEa];
					operandb[m] = operand[m][MEb];
					operanda_[m] = operand[m + corenum][MEa];
					operandb_[m] = operand[m + corenum][MEb];
				}

				for (int m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m],
									 operanda_out[m], operandb_out[m]);

					inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m],
									 operanda_2out[m], operandb_2out[m]);
				}
				j++;

				for (int m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					operand[m][MEa] = static_cast<UDTYPE>(operanda_out[m]) - (modulus & static_cast<UDTYPE>(-static_cast<DTYPE>(operanda_out[m] >= modulus)));
					operand[m][MEb] = static_cast<UDTYPE>(operandb_out[m]) - (modulus & static_cast<UDTYPE>(-static_cast<DTYPE>(operandb_out[m] >= modulus)));
					operand[m + corenum][MEa] = static_cast<UDTYPE>(operanda_2out[m]) - (modulus & static_cast<UDTYPE>(-static_cast<DTYPE>(operanda_2out[m] >= modulus)));
					operand[m + corenum][MEb] = static_cast<UDTYPE>(operandb_2out[m]) - (modulus & static_cast<UDTYPE>(-static_cast<DTYPE>(operandb_2out[m] >= modulus)));
				}

				if (((j + 1) & temp3) == 0)
				{
					MEa += stepsize + 1;
					MEb += stepsize + 1;
				}
				else
				{
					MEa += 1;
					MEb += 1;
				}
			}
		}
	}
}

template <unsigned modcount, unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax, unsigned type>
void intt_4core_mods_new(UDTYPE operand[modcount][bramnum][bramsize],
						 UDTYPE idt[RPBRAMNUM][RPBRAMSIZE],
						 UDTYPE sidt[RPBRAMNUM][RPBRAMSIZE],
						 UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE operanda_out[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_out complete dim = 0
	UDTYPE operandb_out[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_out complete dim = 0
	UDTYPE operanda_2out[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_2out complete dim = 0
	UDTYPE operandb_2out[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_2out complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			// 			for (uint32_t m = 0; m < corenum; m++)
			// 			{
			// #pragma HLS UNROLL
			// 				uint32_t x = (arrayWindex[m] + temp_adder);
			// 				//				w[m]=rp[x];
			// 				//				sw[m]=srp[x];
			// 				w[m] = idt[0];
			// 				sw[m] = sidt[0];
			// 			}

			if (i == 0)
			{
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					IDXTYPE x = (arrayWindex[0] + temp_adder);
					IDXTYPE rp_idx0_j;
					rp_idx0_j = x >> L_RPBRAMNUM;

					w[0] = idt[0][rp_idx0_j];
					sw[0] = sidt[0][rp_idx0_j];
					w[1] = idt[1][rp_idx0_j];
					sw[1] = sidt[1][rp_idx0_j];
					w[2] = idt[2][rp_idx0_j];
					sw[2] = sidt[2][rp_idx0_j];
					w[3] = idt[3][rp_idx0_j];
					sw[3] = sidt[3][rp_idx0_j];

					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][1][MEb];
					operanda[s][1] = operand[s][2][MEa];
					operandb[s][1] = operand[s][3][MEb];
					operanda[s][2] = operand[s][4][MEa];
					operandb[s][2] = operand[s][5][MEb];
					operanda[s][3] = operand[s][6][MEa];
					operandb[s][3] = operand[s][7][MEb];

					for (uint32_t m = 0; m < corenum; m++)
					{
#pragma HLS UNROLL

						inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
										 operanda[s][m], operandb[s][m], operanda_[s][m], operandb_[s][m]);
					}

					operand[s][0][MEa] = operanda_[s][0];
					operand[s][1][MEb] = operandb_[s][0];
					operand[s][2][MEa] = operanda_[s][1];
					operand[s][3][MEb] = operandb_[s][1];
					operand[s][4][MEa] = operanda_[s][2];
					operand[s][5][MEb] = operandb_[s][2];
					operand[s][6][MEa] = operanda_[s][3];
					operand[s][7][MEb] = operandb_[s][3];
				}
			}
			else if (i == 1)
			{
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					IDXTYPE x = (arrayWindex[0] + temp_adder);
					IDXTYPE rp_idx0_i, rp_idx0_j;
					rp_idx0_i = x(L_RPBRAMNUM - 1, 0);
					rp_idx0_j = x >> L_RPBRAMNUM;

					UDTYPE rp_t1, srp_t1;
					UDTYPE rp_t2, srp_t2;

					rp_t1 = idt[rp_idx0_i][rp_idx0_j];
					srp_t1 = sidt[rp_idx0_i][rp_idx0_j];

					rp_t2 = idt[rp_idx0_i + 1][rp_idx0_j];
					srp_t2 = sidt[rp_idx0_i + 1][rp_idx0_j];

					w[0] = rp_t1;
					sw[0] = srp_t1;
					w[1] = rp_t1;
					sw[1] = srp_t1;
					w[2] = rp_t2;
					sw[2] = srp_t2;
					w[3] = rp_t2;
					sw[3] = srp_t2;

					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][2][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][3][MEb];
					operanda[s][2] = operand[s][4][MEa];
					operandb[s][2] = operand[s][6][MEb];
					operanda[s][3] = operand[s][5][MEa];
					operandb[s][3] = operand[s][7][MEb];

					for (uint32_t m = 0; m < corenum; m++)
					{
#pragma HLS UNROLL
						inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
										 operanda[s][m], operandb[s][m], operanda_[s][m], operandb_[s][m]);
					}

					operand[s][0][MEa] = operanda_[s][0];
					operand[s][2][MEb] = operandb_[s][0];
					operand[s][1][MEa] = operanda_[s][1];
					operand[s][3][MEb] = operandb_[s][1];
					operand[s][4][MEa] = operanda_[s][2];
					operand[s][6][MEb] = operandb_[s][2];
					operand[s][5][MEa] = operanda_[s][3];
					operand[s][7][MEb] = operandb_[s][3];
				}
			}
			else if (i == 2)
			{
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					for (uint32_t m = 0; m < corenum; m++)
					{
#pragma HLS UNROLL

						IDXTYPE x = (arrayWindex[0] + temp_adder);
						IDXTYPE rp_idx0_i, rp_idx0_j;
						rp_idx0_i = x(L_RPBRAMNUM - 1, 0);
						rp_idx0_j = x >> L_RPBRAMNUM;
						// cout << "x = " << x << endl;

						w[m] = idt[rp_idx0_i][rp_idx0_j];
						sw[m] = sidt[rp_idx0_i][rp_idx0_j];
					}

					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][4][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][5][MEb];
					operanda[s][2] = operand[s][2][MEa];
					operandb[s][2] = operand[s][6][MEb];
					operanda[s][3] = operand[s][3][MEa];
					operandb[s][3] = operand[s][7][MEb];

					for (uint32_t m = 0; m < corenum; m++)
					{
#pragma HLS UNROLL
						inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
										 operanda[s][m], operandb[s][m], operanda_[s][m], operandb_[s][m]);
					}

					operand[s][0][MEa] = operanda_[s][0];
					operand[s][4][MEb] = operandb_[s][0];
					operand[s][1][MEa] = operanda_[s][1];
					operand[s][5][MEb] = operandb_[s][1];
					operand[s][2][MEa] = operanda_[s][2];
					operand[s][6][MEb] = operandb_[s][2];
					operand[s][3][MEa] = operanda_[s][3];
					operand[s][7][MEb] = operandb_[s][3];
				}
			}

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1_without_last:
	for (uint32_t i = stage_num - stage1_max; i < stage_num - 1; i++)
	{

		uint32_t stepsize = 1 << (i - (stage_num - stage1_max));
		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;
		uint32_t stage_tempa = 12;

		MEa = 0;
		MEb = stepsize;

		// cout << "stepsize" << stepsize << endl;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			// UDTYPE w0 = idt[temp1 + (j >> temp2)];
			// UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			// UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			// UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			IDXTYPE rp_idx0 = temp1 + (j >> temp2);
			IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
			IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

			rp_idx0_i = rp_idx0(L_RPBRAMNUM - 1, 0);
			rp_idx0_j = rp_idx0 >> L_RPBRAMNUM;
			rp_idx00_i = rp_idx00(L_RPBRAMNUM - 1, 0);
			rp_idx00_j = rp_idx00 >> L_RPBRAMNUM;

			UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
			UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
			UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
			UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				for (int m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					operanda[s][m] = operand[s][m][MEa];
					operandb[s][m] = operand[s][m][MEb];
					operanda_[s][m] = operand[s][m + corenum][MEa];
					operandb_[s][m] = operand[s][m + corenum][MEb];
				}
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL

					inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[s][m], operandb[s][m],
									 operanda_out[s][m], operandb_out[s][m]);

					inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[s][m], operandb_[s][m],
									 operanda_2out[s][m], operandb_2out[s][m]);
				}
			}
			j++;

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				for (int m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL

					operand[s][m][MEa] = operanda_out[s][m];
					operand[s][m][MEb] = operandb_out[s][m];
					operand[s][m + corenum][MEa] = operanda_2out[s][m];
					operand[s][m + corenum][MEb] = operandb_2out[s][m];
				}
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}

	if (type == 0)
	{
	type0_intt_stage1_last:
		for (uint32_t i = stage_num - 1; i < stage_num; i++)
		{

			uint32_t stepsize = 1 << (i - (stage_num - stage1_max));
			uint32_t temp2 = i - (stage_num - stage1_max - 1);
			uint32_t temp1 = step_num >> temp2;
			uint32_t temp3 = (one_u64 << temp2) - one_u64;
			uint32_t stage_tempa = 12;

			MEa = 0;
			MEb = stepsize;

			// cout << "stepsize" << stepsize << endl;

		type0_intt_stage1_inner_last:
			for (uint32_t j = 0; j < step_num; j++)
			{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

				// UDTYPE w0 = idt[temp1 + (j >> temp2)];
				// UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
				// UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
				// UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

				IDXTYPE rp_idx0 = temp1 + (j >> temp2);
				IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
				IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

				rp_idx0_i = rp_idx0(L_RPBRAMNUM - 1, 0);
				rp_idx0_j = rp_idx0 >> L_RPBRAMNUM;
				rp_idx00_i = rp_idx00(L_RPBRAMNUM - 1, 0);
				rp_idx00_j = rp_idx00 >> L_RPBRAMNUM;

				UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
				UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
				UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
				UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					for (int m = 0; m < corenum; m++)
					{
#pragma HLS UNROLL
						operanda[s][m] = operand[s][m][MEa];
						operandb[s][m] = operand[s][m][MEb];
						operanda_[s][m] = operand[s][m + corenum][MEa];
						operandb_[s][m] = operand[s][m + corenum][MEb];
					}
				}

				for (int m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					for (int s = 0; s < modcount; s++)
					{
#pragma HLS UNROLL
						inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[s][m], operandb[s][m],
										 operanda_out[s][m], operandb_out[s][m]);

						inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[s][m], operandb_[s][m],
										 operanda_2out[s][m], operandb_2out[s][m]);
					}
				}
				j++;

				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					for (int m = 0; m < corenum; m++)
					{
#pragma HLS UNROLL
						operand[s][m][MEa] = operanda_out[s][m];
						operand[s][m][MEb] = operandb_out[s][m];
						operand[s][m + corenum][MEa] = operanda_2out[s][m];
						operand[s][m + corenum][MEb] = operandb_2out[s][m];
					}
				}

				if (((j + 1) & temp3) == 0)
				{
					MEa += stepsize + 1;
					MEb += stepsize + 1;
				}
				else
				{
					MEa += 1;
					MEb += 1;
				}
			}
		}
	}
	else if (type == 1)
	{
	type1_intt_stage1_last:
		for (uint32_t i = stage_num - 1; i < stage_num; i++)
		{

			uint32_t stepsize = 1 << (i - (stage_num - stage1_max));
			uint32_t temp2 = i - (stage_num - stage1_max - 1);
			uint32_t temp1 = step_num >> temp2;
			uint32_t temp3 = (one_u64 << temp2) - one_u64;
			uint32_t stage_tempa = 12;

			MEa = 0;
			MEb = stepsize;

			// cout << "stepsize" << stepsize << endl;

		type1_intt_stage1_inner_last:
			for (uint32_t j = 0; j < step_num; j++)
			{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

				IDXTYPE rp_idx0 = temp1 + (j >> temp2);
				IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
				IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

				rp_idx0_i = rp_idx0(L_RPBRAMNUM - 1, 0);
				rp_idx0_j = rp_idx0 >> L_RPBRAMNUM;
				rp_idx00_i = rp_idx00(L_RPBRAMNUM - 1, 0);
				rp_idx00_j = rp_idx00 >> L_RPBRAMNUM;

				UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
				UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
				UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
				UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					for (int m = 0; m < corenum; m++)
					{
#pragma HLS UNROLL
						operanda[s][m] = operand[s][m][MEa];
						operandb[s][m] = operand[s][m][MEb];
						operanda_[s][m] = operand[s][m + corenum][MEa];
						operandb_[s][m] = operand[s][m + corenum][MEb];
					}
				}

				for (int m = 0; m < corenum; m++)
				{
#pragma HLS UNROLL
					for (int s = 0; s < modcount; s++)
					{
#pragma HLS UNROLL
						inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[s][m], operandb[s][m],
										 operanda_out[s][m], operandb_out[s][m]);

						inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[s][m], operandb_[s][m],
										 operanda_2out[s][m], operandb_2out[s][m]);
					}
				}
				j++;

				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					for (int m = 0; m < corenum; m++)
					{
#pragma HLS UNROLL
						operand[s][m][MEa] = static_cast<UDTYPE>(operanda_out[s][m]) - (modulus & static_cast<UDTYPE>(-static_cast<DTYPE>(operanda_out[s][m] >= modulus)));
						operand[s][m][MEb] = static_cast<UDTYPE>(operandb_out[s][m]) - (modulus & static_cast<UDTYPE>(-static_cast<DTYPE>(operandb_out[s][m] >= modulus)));
						operand[s][m + corenum][MEa] = static_cast<UDTYPE>(operanda_2out[s][m]) - (modulus & static_cast<UDTYPE>(-static_cast<DTYPE>(operanda_2out[s][m] >= modulus)));
						operand[s][m + corenum][MEb] = static_cast<UDTYPE>(operandb_2out[s][m]) - (modulus & static_cast<UDTYPE>(-static_cast<DTYPE>(operandb_2out[s][m] >= modulus)));
					}
				}

				if (((j + 1) & temp3) == 0)
				{
					MEa += stepsize + 1;
					MEb += stepsize + 1;
				}
				else
				{
					MEa += 1;
					MEb += 1;
				}
			}
		}
	}
}

template <unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax>
void intt_8core(UDTYPE operand[bramnum][bramsize],
				UDTYPE idt[N], UDTYPE sidt[N], UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				uint32_t x = (arrayWindex[m] + temp_adder);
				//				w[m]=rp[x];
				//				sw[m]=srp[x];
				w[m] = idt[0];
				sw[m] = sidt[0];
			}

			if (i == 0)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[1][MEb];
				operanda[1] = operand[2][MEa];
				operandb[1] = operand[3][MEb];
				operanda[2] = operand[4][MEa];
				operandb[2] = operand[5][MEb];
				operanda[3] = operand[6][MEa];
				operandb[3] = operand[7][MEb];
				operanda[4] = operand[8][MEa];
				operandb[4] = operand[9][MEb];
				operanda[5] = operand[10][MEa];
				operandb[5] = operand[11][MEb];
				operanda[6] = operand[12][MEa];
				operandb[6] = operand[13][MEb];
				operanda[7] = operand[14][MEa];
				operandb[7] = operand[15][MEb];
			}
			else if (i == 1)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[2][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[3][MEb];
				operanda[2] = operand[4][MEa];
				operandb[2] = operand[6][MEb];
				operanda[3] = operand[5][MEa];
				operandb[3] = operand[7][MEb];
				operanda[4] = operand[8][MEa];
				operandb[4] = operand[10][MEb];
				operanda[5] = operand[9][MEa];
				operandb[5] = operand[11][MEb];
				operanda[6] = operand[12][MEa];
				operandb[6] = operand[14][MEb];
				operanda[7] = operand[13][MEa];
				operandb[7] = operand[15][MEb];
			}
			else if (i == 2)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[4][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[5][MEb];
				operanda[2] = operand[2][MEa];
				operandb[2] = operand[6][MEb];
				operanda[3] = operand[3][MEa];
				operandb[3] = operand[7][MEb];
				operanda[4] = operand[8][MEa];
				operandb[4] = operand[12][MEb];
				operanda[5] = operand[9][MEa];
				operandb[5] = operand[13][MEb];
				operanda[6] = operand[10][MEa];
				operandb[6] = operand[14][MEb];
				operanda[7] = operand[11][MEa];
				operandb[7] = operand[15][MEb];
			}
			else if (i == 3)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[8][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[9][MEb];
				operanda[2] = operand[2][MEa];
				operandb[2] = operand[10][MEb];
				operanda[3] = operand[3][MEa];
				operandb[3] = operand[11][MEb];
				operanda[4] = operand[4][MEa];
				operandb[4] = operand[12][MEb];
				operanda[5] = operand[5][MEa];
				operandb[5] = operand[13][MEb];
				operanda[6] = operand[6][MEa];
				operandb[6] = operand[14][MEb];
				operanda[7] = operand[7][MEa];
				operandb[7] = operand[15][MEb];
			}

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w[m], sw[m], modulus, two_times_modulus, operanda[m], operandb[m]);
			}

			if (i == 0)
			{
				operand[0][MEa] = operanda[0];
				operand[1][MEb] = operandb[0];
				operand[2][MEa] = operanda[1];
				operand[3][MEb] = operandb[1];
				operand[4][MEa] = operanda[2];
				operand[5][MEb] = operandb[2];
				operand[6][MEa] = operanda[3];
				operand[7][MEb] = operandb[3];
				operand[8][MEa] = operanda[4];
				operand[9][MEb] = operandb[4];
				operand[10][MEa] = operanda[5];
				operand[11][MEb] = operandb[5];
				operand[12][MEa] = operanda[6];
				operand[13][MEb] = operandb[6];
				operand[14][MEa] = operanda[7];
				operand[15][MEb] = operandb[7];
			}
			else if (i == 1)
			{
				operand[0][MEa] = operanda[0];
				operand[2][MEb] = operandb[0];
				operand[1][MEa] = operanda[1];
				operand[3][MEb] = operandb[1];
				operand[4][MEa] = operanda[2];
				operand[6][MEb] = operandb[2];
				operand[5][MEa] = operanda[3];
				operand[7][MEb] = operandb[3];
				operand[8][MEa] = operanda[4];
				operand[10][MEb] = operandb[4];
				operand[9][MEa] = operanda[5];
				operand[11][MEb] = operandb[5];
				operand[12][MEa] = operanda[6];
				operand[14][MEb] = operandb[6];
				operand[13][MEa] = operanda[7];
				operand[15][MEb] = operandb[7];
			}
			else if (i == 2)
			{
				operand[0][MEa] = operanda[0];
				operand[4][MEb] = operandb[0];
				operand[1][MEa] = operanda[1];
				operand[5][MEb] = operandb[1];
				operand[2][MEa] = operanda[2];
				operand[6][MEb] = operandb[2];
				operand[3][MEa] = operanda[3];
				operand[7][MEb] = operandb[3];
				operand[8][MEa] = operanda[4];
				operand[12][MEb] = operandb[4];
				operand[9][MEa] = operanda[5];
				operand[13][MEb] = operandb[5];
				operand[10][MEa] = operanda[6];
				operand[14][MEb] = operandb[6];
				operand[11][MEa] = operanda[7];
				operand[15][MEb] = operandb[7];
			}
			else if (i == 3)
			{
				operand[0][MEa] = operanda[0];
				operand[8][MEb] = operandb[0];
				operand[1][MEa] = operanda[1];
				operand[9][MEb] = operandb[1];
				operand[2][MEa] = operanda[2];
				operand[10][MEb] = operandb[2];
				operand[3][MEa] = operanda[3];
				operand[11][MEb] = operandb[3];
				operand[4][MEa] = operanda[4];
				operand[12][MEb] = operandb[4];
				operand[5][MEa] = operanda[5];
				operand[13][MEb] = operandb[5];
				operand[6][MEa] = operanda[6];
				operand[14][MEb] = operandb[6];
				operand[7][MEa] = operanda[7];
				operand[15][MEb] = operandb[7];
			}

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1:
	for (uint32_t i = stage_num - stage1_max; i < stage_num; i++)
	{
		uint32_t stepsize = step_num >> (i + 1);

		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			UDTYPE w0 = idt[temp1 + (j >> temp2)];
			UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operanda[m] = operand[m][MEa];
				operandb[m] = operand[m][MEb];
				operanda_[m] = operand[m + corenum][MEa];
				operandb_[m] = operand[m + corenum][MEb];
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m]);
			}
			j++;
			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m]);
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operand[m][MEa] = operanda[m];
				operand[m][MEb] = operandb[m];
				operand[m + corenum][MEa] = operanda_[m];
				operand[m + corenum][MEb] = operandb_[m];
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}
}

template <unsigned total_poly, unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax, unsigned type>
void intt_8core_new(UDTYPE operand[bramnum][bramsize],
					UDTYPE idt[RPBRAMNUM][RPBRAMSIZE],
					UDTYPE sidt[RPBRAMNUM][RPBRAMSIZE],
					UDTYPE modulus)
{
#pragma HLS INLINE off

    UDTYPE operanda[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
    UDTYPE operandb[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
    UDTYPE operanda_[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
    UDTYPE operandb_[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

    UDTYPE operanda_out[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_out complete dim = 0
    UDTYPE operandb_out[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_out complete dim = 0
    UDTYPE operanda_2out[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_2out complete dim = 0
    UDTYPE operandb_2out[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_2out complete dim = 0

    UDTYPE two_times_modulus = modulus << 1;
    UDTYPE one_u64 = 1;
    uint32_t step_num = bramsize;
    uint32_t stage_num = STAGENUM;
    uint32_t stage1_max = stagemax;

    uint32_t MEa;
    uint32_t MEb;

    uint32_t arrayWindex;
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

    intt_stage2:
    for (uint32_t i = 0; i < stage_num - stage1_max; i++) {

        uint32_t temp4 = stage_num - 1 - i;
        uint32_t temp_par = corenum >> i;
        // uint32_t temp5 = i - stage1_max;

        MEa = 0;
        MEb = 0;

        arrayWindex = (1 << temp4);

        intt_stage2_inner:
        for (uint32_t j = 0; j < step_num; j++) {
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

            UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
            UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

            uint32_t temp_adder = temp_par * j;


            if (i == 0) {
                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    IDXTYPE x = (arrayWindex + temp_adder);
                    IDXTYPE rp_idx0_j;
                    rp_idx0_j = x >> logRPBRAMNUM;

                    w[0] = idt[0][rp_idx0_j];
                    sw[0] = sidt[0][rp_idx0_j];
                    w[1] = idt[1][rp_idx0_j];
                    sw[1] = sidt[1][rp_idx0_j];
                    w[2] = idt[2][rp_idx0_j];
                    sw[2] = sidt[2][rp_idx0_j];
                    w[3] = idt[3][rp_idx0_j];
                    sw[3] = sidt[3][rp_idx0_j];

                    w[4] = idt[0][rp_idx0_j];
                    w[5] = idt[1][rp_idx0_j];
                    w[6] = idt[2][rp_idx0_j];
                    w[7] = idt[3][rp_idx0_j];
                    sw[4] = sidt[0][rp_idx0_j];
                    sw[5] = sidt[1][rp_idx0_j];
                    sw[6] = sidt[2][rp_idx0_j];
                    sw[7] = sidt[3][rp_idx0_j];

                    operanda[0] = operand[0][MEa];
                    operandb[0] = operand[1][MEb];
                    operanda[1] = operand[2][MEa];
                    operandb[1] = operand[3][MEb];
                    operanda[2] = operand[4][MEa];
                    operandb[2] = operand[5][MEb];
                    operanda[3] = operand[6][MEa];
                    operandb[3] = operand[7][MEb];
                    operanda[4] = operand[8][MEa];
                    operandb[4] = operand[9][MEb];
                    operanda[5] = operand[10][MEa];
                    operandb[5] = operand[11][MEb];
                    operanda[6] = operand[12][MEa];
                    operandb[6] = operand[13][MEb];
                    operanda[7] = operand[14][MEa];
                    operandb[7] = operand[15][MEb];

                    for (uint32_t m = 0; m < corenum; m++) {
#pragma HLS UNROLL

                        inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
                                         operanda[m], operandb[m], operanda_[m], operandb_[m]);
                    }

                    operand[0][MEa] = operanda_[0];
                    operand[1][MEb] = operandb_[0];
                    operand[2][MEa] = operanda_[1];
                    operand[3][MEb] = operandb_[1];
                    operand[4][MEa] = operanda_[2];
                    operand[5][MEb] = operandb_[2];
                    operand[6][MEa] = operanda_[3];
                    operand[7][MEb] = operandb_[3];
                    operand[8][MEa] = operanda_[4];
                    operand[9][MEb] = operandb_[4];
                    operand[10][MEa] = operanda_[5];
                    operand[11][MEb] = operandb_[5];
                    operand[12][MEa] = operanda_[6];
                    operand[13][MEb] = operandb_[6];
                    operand[14][MEa] = operanda_[7];
                    operand[15][MEb] = operandb_[7];
                }
            }
            else if (i == 1) {
                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    IDXTYPE x = (arrayWindex + temp_adder);
                    IDXTYPE rp_idx0_i, rp_idx0_j;
                    rp_idx0_i = x(logRPBRAMNUM - 1, 0);
                    rp_idx0_j = x >> logRPBRAMNUM;

                    UDTYPE rp_t1, srp_t1;
                    UDTYPE rp_t2, srp_t2;

                    rp_t1 = idt[rp_idx0_i][rp_idx0_j];
                    srp_t1 = sidt[rp_idx0_i][rp_idx0_j];

                    rp_t2 = idt[rp_idx0_i + 1][rp_idx0_j];
                    srp_t2 = sidt[rp_idx0_i + 1][rp_idx0_j];

                    w[0] = rp_t1;
                    sw[0] = srp_t1;
                    w[1] = rp_t1;
                    sw[1] = srp_t1;
                    w[2] = rp_t1;
                    sw[2] = srp_t1;
                    w[3] = rp_t1;
                    sw[3] = srp_t1;

                    w[4] = rp_t2;
                    w[5] = rp_t2;
                    w[6] = rp_t2;
                    w[7] = rp_t2;
                    sw[4] = srp_t2;
                    sw[5] = srp_t2;
                    sw[6] = srp_t2;
                    sw[7] = srp_t2;

                    operanda[0] = operand[0][MEa];
                    operandb[0] = operand[2][MEb];
                    operanda[1] = operand[1][MEa];
                    operandb[1] = operand[3][MEb];
                    operanda[2] = operand[4][MEa];
                    operandb[2] = operand[6][MEb];
                    operanda[3] = operand[5][MEa];
                    operandb[3] = operand[7][MEb];
                    operanda[4] = operand[8][MEa];
                    operandb[4] = operand[10][MEb];
                    operanda[5] = operand[9][MEa];
                    operandb[5] = operand[11][MEb];
                    operanda[6] = operand[12][MEa];
                    operandb[6] = operand[14][MEb];
                    operanda[7] = operand[13][MEa];
                    operandb[7] = operand[15][MEb];

                    for (uint32_t m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
                                         operanda[m], operandb[m], operanda_[m], operandb_[m]);
                    }

                    operand[0][MEa] = operanda_[0];
                    operand[2][MEb] = operandb_[0];
                    operand[1][MEa] = operanda_[1];
                    operand[3][MEb] = operandb_[1];
                    operand[4][MEa] = operanda_[2];
                    operand[6][MEb] = operandb_[2];
                    operand[5][MEa] = operanda_[3];
                    operand[7][MEb] = operandb_[3];
                    operand[8][MEa] = operanda_[4];
                    operand[10][MEb] = operandb_[4];
                    operand[9][MEa] = operanda_[5];
                    operand[11][MEb] = operandb_[5];
                    operand[12][MEa] = operanda_[6];
                    operand[14][MEb] = operandb_[6];
                    operand[13][MEa] = operanda_[7];
                    operand[15][MEb] = operandb_[7];
                }
            }
            else if (i == 2) {
                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL

                    IDXTYPE x = (arrayWindex + temp_adder);
                    IDXTYPE rp_idx0_i, rp_idx0_j;
                    rp_idx0_i = x(logRPBRAMNUM - 1, 0);
                    rp_idx0_j = x >> logRPBRAMNUM;
                    // cout << "x = " << x << endl;

                    UDTYPE rp_t1, srp_t1;
                    UDTYPE rp_t2, srp_t2;
                    UDTYPE rp_t3, srp_t3;
                    UDTYPE rp_t4, srp_t4;

                    rp_t1 = idt[rp_idx0_i][rp_idx0_j];
                    srp_t1 = sidt[rp_idx0_i][rp_idx0_j];

                    rp_t2 = idt[rp_idx0_i + 1][rp_idx0_j];
                    srp_t2 = sidt[rp_idx0_i + 1][rp_idx0_j];

                    rp_t3 = idt[rp_idx0_i + 2][rp_idx0_j];
                    srp_t3 = sidt[rp_idx0_i + 2][rp_idx0_j];

                    rp_t4 = idt[rp_idx0_i + 3][rp_idx0_j];
                    srp_t4 = sidt[rp_idx0_i + 3][rp_idx0_j];

                    w[0] = rp_t1;
                    sw[0] = srp_t1;
                    w[1] = rp_t1;
                    sw[1] = srp_t1;

                    w[2] = rp_t2;
                    sw[2] = srp_t2;
                    w[3] = rp_t2;
                    sw[3] = srp_t2;

                    w[4] = rp_t3;
                    w[5] = rp_t3;
                    w[6] = rp_t4;
                    w[7] = rp_t4;
                    sw[4] = srp_t3;
                    sw[5] = srp_t3;
                    sw[6] = srp_t4;
                    sw[7] = srp_t4;

                    operanda[0] = operand[0][MEa];
                    operandb[0] = operand[4][MEb];
                    operanda[1] = operand[1][MEa];
                    operandb[1] = operand[5][MEb];
                    operanda[2] = operand[2][MEa];
                    operandb[2] = operand[6][MEb];
                    operanda[3] = operand[3][MEa];
                    operandb[3] = operand[7][MEb];
                    operanda[4] = operand[8][MEa];
                    operandb[4] = operand[12][MEb];
                    operanda[5] = operand[9][MEa];
                    operandb[5] = operand[13][MEb];
                    operanda[6] = operand[10][MEa];
                    operandb[6] = operand[14][MEb];
                    operanda[7] = operand[11][MEa];
                    operandb[7] = operand[15][MEb];

                    for (uint32_t m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
                                         operanda[m], operandb[m], operanda_[m], operandb_[m]);
                    }

                    operand[0][MEa] = operanda_[0];
                    operand[4][MEb] = operandb_[0];
                    operand[1][MEa] = operanda_[1];
                    operand[5][MEb] = operandb_[1];
                    operand[2][MEa] = operanda_[2];
                    operand[6][MEb] = operandb_[2];
                    operand[3][MEa] = operanda_[3];
                    operand[7][MEb] = operandb_[3];
                    operand[8][MEa] = operanda_[4];
                    operand[12][MEb] = operandb_[4];
                    operand[9][MEa] = operanda_[5];
                    operand[13][MEb] = operandb_[5];
                    operand[10][MEa] = operanda_[6];
                    operand[14][MEb] = operandb_[6];
                    operand[11][MEa] = operanda_[7];
                    operand[15][MEb] = operandb_[7];
                }
            }
            else if (i == 3) {
                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL

                    for (uint32_t m = 0; m < corenum; m++) {
#pragma HLS UNROLL

                        IDXTYPE x = (arrayWindex + temp_adder);
                        IDXTYPE rp_idx0_i, rp_idx0_j;
                        rp_idx0_i = x(logRPBRAMNUM - 1, 0);
                        rp_idx0_j = x >> logRPBRAMNUM;
                        // cout << "x = " << x << endl;

                        w[m] = idt[rp_idx0_i][rp_idx0_j];
                        sw[m] = sidt[rp_idx0_i][rp_idx0_j];
                    }

                    operanda[0] = operand[0][MEa];
                    operandb[0] = operand[8][MEb];
                    operanda[1] = operand[1][MEa];
                    operandb[1] = operand[9][MEb];
                    operanda[2] = operand[2][MEa];
                    operandb[2] = operand[10][MEb];
                    operanda[3] = operand[3][MEa];
                    operandb[3] = operand[11][MEb];
                    operanda[4] = operand[4][MEa];
                    operandb[4] = operand[12][MEb];
                    operanda[5] = operand[5][MEa];
                    operandb[5] = operand[13][MEb];
                    operanda[6] = operand[6][MEa];
                    operandb[6] = operand[14][MEb];
                    operanda[7] = operand[7][MEa];
                    operandb[7] = operand[15][MEb];

                    for (uint32_t m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
                                         operanda[m], operandb[m], operanda_[m], operandb_[m]);
                    }

                    operand[0][MEa] = operanda_[0];
                    operand[8][MEb] = operandb_[0];
                    operand[1][MEa] = operanda_[1];
                    operand[9][MEb] = operandb_[1];
                    operand[2][MEa] = operanda_[2];
                    operand[10][MEb] = operandb_[2];
                    operand[3][MEa] = operanda_[3];
                    operand[11][MEb] = operandb_[3];
                    operand[4][MEa] = operanda_[4];
                    operand[12][MEb] = operandb_[4];
                    operand[5][MEa] = operanda_[5];
                    operand[13][MEb] = operandb_[5];
                    operand[6][MEa] = operanda_[6];
                    operand[14][MEb] = operandb_[6];
                    operand[7][MEa] = operanda_[7];
                    operand[15][MEb] = operandb_[7];
                }
            }

            MEa += 1;
            MEb += 1;
        }
    }

    intt_stage1_without_last:
    for (uint32_t i = stage_num - stage1_max; i < stage_num - 1; i++) {

        uint32_t stepsize = 1 << (i - (stage_num - stage1_max));
        uint32_t temp2 = i - (stage_num - stage1_max - 1);
        uint32_t temp1 = step_num >> temp2;
        uint32_t temp3 = (one_u64 << temp2) - one_u64;
        uint32_t stage_tempa = 12;

        MEa = 0;
        MEb = stepsize;

        // cout << "stepsize" << stepsize << endl;

        intt_stage1_inner:
        for (uint32_t j = 0; j < step_num; j++) {
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

            // UDTYPE w0 = idt[temp1 + (j >> temp2)];
            // UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
            // UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
            // UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

            IDXTYPE rp_idx0 = temp1 + (j >> temp2);
            IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
            IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

            rp_idx0_i = rp_idx0(logRPBRAMNUM - 1, 0);
            rp_idx0_j = rp_idx0 >> logRPBRAMNUM;
            rp_idx00_i = rp_idx00(logRPBRAMNUM - 1, 0);
            rp_idx00_j = rp_idx00 >> logRPBRAMNUM;

            UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
            UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
            UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
            UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

            for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                    operanda[m] = operand[m][MEa];
                    operandb[m] = operand[m][MEb];
                    operanda_[m] = operand[m + corenum][MEa];
                    operandb_[m] = operand[m + corenum][MEb];
                }
            }

            for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL

                    inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m],
                                     operanda_out[m], operandb_out[m]);

                    inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m],
                                     operanda_2out[m], operandb_2out[m]);
                }
            }
            j++;

            for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL

                    operand[m][MEa] = operanda_out[m];
                    operand[m][MEb] = operandb_out[m];
                    operand[m + corenum][MEa] = operanda_2out[m];
                    operand[m + corenum][MEb] = operandb_2out[m];
                }
            }

            if (((j + 1) & temp3) == 0) {
                MEa += stepsize + 1;
                MEb += stepsize + 1;
            }
            else {
                MEa += 1;
                MEb += 1;
            }
        }
    }

    if (type == 0) {
        type0_intt_stage1_last:
        for (uint32_t i = stage_num - 1; i < stage_num; i++) {

            uint32_t stepsize = 1 << (i - (stage_num - stage1_max));
            uint32_t temp2 = i - (stage_num - stage1_max - 1);
            uint32_t temp1 = step_num >> temp2;
            uint32_t temp3 = (one_u64 << temp2) - one_u64;
            uint32_t stage_tempa = 12;

            MEa = 0;
            MEb = stepsize;

            // cout << "stepsize" << stepsize << endl;

            type0_intt_stage1_inner_last:
            for (uint32_t j = 0; j < step_num; j++) {
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

                // UDTYPE w0 = idt[temp1 + (j >> temp2)];
                // UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
                // UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
                // UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

                IDXTYPE rp_idx0 = temp1 + (j >> temp2);
                IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
                IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

                rp_idx0_i = rp_idx0(logRPBRAMNUM - 1, 0);
                rp_idx0_j = rp_idx0 >> logRPBRAMNUM;
                rp_idx00_i = rp_idx00(logRPBRAMNUM - 1, 0);
                rp_idx00_j = rp_idx00 >> logRPBRAMNUM;

                UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
                UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
                UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
                UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        operanda[m] = operand[m][MEa];
                        operandb[m] = operand[m][MEb];
                        operanda_[m] = operand[m + corenum][MEa];
                        operandb_[m] = operand[m + corenum][MEb];
                    }
                }

                for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                    for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                        inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m],
                                         operanda_out[m], operandb_out[m]);

                        inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m],
                                         operanda_2out[m], operandb_2out[m]);
                    }
                }
                j++;

                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        operand[m][MEa] = operanda_out[m];
                        operand[m][MEb] = operandb_out[m];
                        operand[m + corenum][MEa] = operanda_2out[m];
                        operand[m + corenum][MEb] = operandb_2out[m];
                    }
                }

                if (((j + 1) & temp3) == 0) {
                    MEa += stepsize + 1;
                    MEb += stepsize + 1;
                }
                else {
                    MEa += 1;
                    MEb += 1;
                }
            }
        }
    }
    else if (type == 1) {
        type1_intt_stage1_last:
        for (uint32_t i = stage_num - 1; i < stage_num; i++) {

            uint32_t stepsize = 1 << (i - (stage_num - stage1_max));
            uint32_t temp2 = i - (stage_num - stage1_max - 1);
            uint32_t temp1 = step_num >> temp2;
            uint32_t temp3 = (one_u64 << temp2) - one_u64;
            uint32_t stage_tempa = 12;

            MEa = 0;
            MEb = stepsize;

            // cout << "stepsize" << stepsize << endl;

            type1_intt_stage1_inner_last:
            for (uint32_t j = 0; j < step_num; j++) {
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

                IDXTYPE rp_idx0 = temp1 + (j >> temp2);
                IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
                IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

                rp_idx0_i = rp_idx0(logRPBRAMNUM - 1, 0);
                rp_idx0_j = rp_idx0 >> logRPBRAMNUM;
                rp_idx00_i = rp_idx00(logRPBRAMNUM - 1, 0);
                rp_idx00_j = rp_idx00 >> logRPBRAMNUM;

                UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
                UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
                UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
                UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        operanda[m] = operand[m][MEa];
                        operandb[m] = operand[m][MEb];
                        operanda_[m] = operand[m + corenum][MEa];
                        operandb_[m] = operand[m + corenum][MEb];
                    }
                }

                for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                    for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                        inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m],
                                         operanda_out[m], operandb_out[m]);

                        inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m],
                                         operanda_2out[m], operandb_2out[m]);
                    }
                }
                j++;

                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        operand[m][MEa] = static_cast<UDTYPE>(operanda_out[m]) - (modulus &
                                                                                  static_cast<UDTYPE>(-static_cast<DTYPE>(
                                                                                          operanda_out[m] >= modulus)));
                        operand[m][MEb] = static_cast<UDTYPE>(operandb_out[m]) - (modulus &
                                                                                  static_cast<UDTYPE>(-static_cast<DTYPE>(
                                                                                          operandb_out[m] >= modulus)));
                        operand[m + corenum][MEa] = static_cast<UDTYPE>(operanda_2out[m]) - (modulus &
                                                                                             static_cast<UDTYPE>(-static_cast<DTYPE>(
                                                                                                     operanda_2out[m] >=
                                                                                                     modulus)));
                        operand[m + corenum][MEb] = static_cast<UDTYPE>(operandb_2out[m]) - (modulus &
                                                                                             static_cast<UDTYPE>(-static_cast<DTYPE>(
                                                                                                     operandb_2out[m] >=
                                                                                                     modulus)));
                    }
                }

                if (((j + 1) & temp3) == 0) {
                    MEa += stepsize + 1;
                    MEb += stepsize + 1;
                }
                else {
                    MEa += 1;
                    MEb += 1;
                }
            }
        }
    }
}



template <unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax>
void intt_16core(UDTYPE operand[bramnum][bramsize],
				 UDTYPE idt[N], UDTYPE sidt[N], UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				uint32_t x = (arrayWindex[m] + temp_adder);
				//				w[m]=rp[x];
				//				sw[m]=srp[x];
				w[m] = idt[0];
				sw[m] = sidt[0];
			}

			if (i == 0)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[1][MEb];
				operanda[1] = operand[2][MEa];
				operandb[1] = operand[3][MEb];
				operanda[2] = operand[4][MEa];
				operandb[2] = operand[5][MEb];
				operanda[3] = operand[6][MEa];
				operandb[3] = operand[7][MEb];
				operanda[4] = operand[8][MEa];
				operandb[4] = operand[9][MEb];
				operanda[5] = operand[10][MEa];
				operandb[5] = operand[11][MEb];
				operanda[6] = operand[12][MEa];
				operandb[6] = operand[13][MEb];
				operanda[7] = operand[14][MEa];
				operandb[7] = operand[15][MEb];
				operanda[8] = operand[16][MEa];
				operandb[8] = operand[17][MEb];
				operanda[9] = operand[18][MEa];
				operandb[9] = operand[19][MEb];
				operanda[10] = operand[20][MEa];
				operandb[10] = operand[21][MEb];
				operanda[11] = operand[22][MEa];
				operandb[11] = operand[23][MEb];
				operanda[12] = operand[24][MEa];
				operandb[12] = operand[25][MEb];
				operanda[13] = operand[26][MEa];
				operandb[13] = operand[27][MEb];
				operanda[14] = operand[28][MEa];
				operandb[14] = operand[29][MEb];
				operanda[15] = operand[30][MEa];
				operandb[15] = operand[31][MEb];
			}
			else if (i == 1)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[2][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[3][MEb];
				operanda[2] = operand[4][MEa];
				operandb[2] = operand[6][MEb];
				operanda[3] = operand[5][MEa];
				operandb[3] = operand[7][MEb];
				operanda[4] = operand[8][MEa];
				operandb[4] = operand[10][MEb];
				operanda[5] = operand[9][MEa];
				operandb[5] = operand[11][MEb];
				operanda[6] = operand[12][MEa];
				operandb[6] = operand[14][MEb];
				operanda[7] = operand[13][MEa];
				operandb[7] = operand[15][MEb];
				operanda[8] = operand[16][MEa];
				operandb[8] = operand[18][MEb];
				operanda[9] = operand[17][MEa];
				operandb[9] = operand[19][MEb];
				operanda[10] = operand[20][MEa];
				operandb[10] = operand[22][MEb];
				operanda[11] = operand[21][MEa];
				operandb[11] = operand[23][MEb];
				operanda[12] = operand[24][MEa];
				operandb[12] = operand[26][MEb];
				operanda[13] = operand[25][MEa];
				operandb[13] = operand[27][MEb];
				operanda[14] = operand[28][MEa];
				operandb[14] = operand[30][MEb];
				operanda[15] = operand[29][MEa];
				operandb[15] = operand[31][MEb];
			}
			else if (i == 2)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[4][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[5][MEb];
				operanda[2] = operand[2][MEa];
				operandb[2] = operand[6][MEb];
				operanda[3] = operand[3][MEa];
				operandb[3] = operand[7][MEb];
				operanda[4] = operand[8][MEa];
				operandb[4] = operand[12][MEb];
				operanda[5] = operand[9][MEa];
				operandb[5] = operand[13][MEb];
				operanda[6] = operand[10][MEa];
				operandb[6] = operand[14][MEb];
				operanda[7] = operand[11][MEa];
				operandb[7] = operand[15][MEb];
				operanda[8] = operand[16][MEa];
				operandb[8] = operand[20][MEb];
				operanda[9] = operand[17][MEa];
				operandb[9] = operand[21][MEb];
				operanda[10] = operand[18][MEa];
				operandb[10] = operand[22][MEb];
				operanda[11] = operand[19][MEa];
				operandb[11] = operand[23][MEb];
				operanda[12] = operand[24][MEa];
				operandb[12] = operand[28][MEb];
				operanda[13] = operand[25][MEa];
				operandb[13] = operand[29][MEb];
				operanda[14] = operand[26][MEa];
				operandb[14] = operand[30][MEb];
				operanda[15] = operand[27][MEa];
				operandb[15] = operand[31][MEb];
			}
			else if (i == 3)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[8][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[9][MEb];
				operanda[2] = operand[2][MEa];
				operandb[2] = operand[10][MEb];
				operanda[3] = operand[3][MEa];
				operandb[3] = operand[11][MEb];
				operanda[4] = operand[4][MEa];
				operandb[4] = operand[12][MEb];
				operanda[5] = operand[5][MEa];
				operandb[5] = operand[13][MEb];
				operanda[6] = operand[6][MEa];
				operandb[6] = operand[14][MEb];
				operanda[7] = operand[7][MEa];
				operandb[7] = operand[15][MEb];
				operanda[8] = operand[16][MEa];
				operandb[8] = operand[24][MEb];
				operanda[9] = operand[17][MEa];
				operandb[9] = operand[25][MEb];
				operanda[10] = operand[18][MEa];
				operandb[10] = operand[26][MEb];
				operanda[11] = operand[19][MEa];
				operandb[11] = operand[27][MEb];
				operanda[12] = operand[20][MEa];
				operandb[12] = operand[28][MEb];
				operanda[13] = operand[21][MEa];
				operandb[13] = operand[29][MEb];
				operanda[14] = operand[22][MEa];
				operandb[14] = operand[30][MEb];
				operanda[15] = operand[23][MEa];
				operandb[15] = operand[31][MEb];
			}
			else if (i == 4)
			{
				operanda[0] = operand[0][MEa];
				operandb[0] = operand[16][MEb];
				operanda[1] = operand[1][MEa];
				operandb[1] = operand[17][MEb];
				operanda[2] = operand[2][MEa];
				operandb[2] = operand[18][MEb];
				operanda[3] = operand[3][MEa];
				operandb[3] = operand[19][MEb];
				operanda[4] = operand[4][MEa];
				operandb[4] = operand[20][MEb];
				operanda[5] = operand[5][MEa];
				operandb[5] = operand[21][MEb];
				operanda[6] = operand[6][MEa];
				operandb[6] = operand[22][MEb];
				operanda[7] = operand[7][MEa];
				operandb[7] = operand[23][MEb];
				operanda[8] = operand[8][MEa];
				operandb[8] = operand[24][MEb];
				operanda[9] = operand[9][MEa];
				operandb[9] = operand[25][MEb];
				operanda[10] = operand[10][MEa];
				operandb[10] = operand[26][MEb];
				operanda[11] = operand[11][MEa];
				operandb[11] = operand[27][MEb];
				operanda[12] = operand[12][MEa];
				operandb[12] = operand[28][MEb];
				operanda[13] = operand[13][MEa];
				operandb[13] = operand[29][MEb];
				operanda[14] = operand[14][MEa];
				operandb[14] = operand[30][MEb];
				operanda[15] = operand[15][MEa];
				operandb[15] = operand[31][MEb];
			}

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w[m], sw[m], modulus, two_times_modulus, operanda[m], operandb[m]);
			}

			if (i == 0)
			{
				operand[0][MEa] = operanda[0];
				operand[1][MEb] = operandb[0];
				operand[2][MEa] = operanda[1];
				operand[3][MEb] = operandb[1];
				operand[4][MEa] = operanda[2];
				operand[5][MEb] = operandb[2];
				operand[6][MEa] = operanda[3];
				operand[7][MEb] = operandb[3];
				operand[8][MEa] = operanda[4];
				operand[9][MEb] = operandb[4];
				operand[10][MEa] = operanda[5];
				operand[11][MEb] = operandb[5];
				operand[12][MEa] = operanda[6];
				operand[13][MEb] = operandb[6];
				operand[14][MEa] = operanda[7];
				operand[15][MEb] = operandb[7];
				operand[16][MEa] = operanda[8];
				operand[17][MEb] = operandb[8];
				operand[18][MEa] = operanda[9];
				operand[19][MEb] = operandb[9];
				operand[20][MEa] = operanda[10];
				operand[21][MEb] = operandb[10];
				operand[22][MEa] = operanda[11];
				operand[23][MEb] = operandb[11];
				operand[24][MEa] = operanda[12];
				operand[25][MEb] = operandb[12];
				operand[26][MEa] = operanda[13];
				operand[27][MEb] = operandb[13];
				operand[28][MEa] = operanda[14];
				operand[29][MEb] = operandb[14];
				operand[30][MEa] = operanda[15];
				operand[31][MEb] = operandb[15];
			}
			else if (i == 1)
			{
				operand[0][MEa] = operanda[0];
				operand[2][MEb] = operandb[0];
				operand[1][MEa] = operanda[1];
				operand[3][MEb] = operandb[1];
				operand[4][MEa] = operanda[2];
				operand[6][MEb] = operandb[2];
				operand[5][MEa] = operanda[3];
				operand[7][MEb] = operandb[3];
				operand[8][MEa] = operanda[4];
				operand[10][MEb] = operandb[4];
				operand[9][MEa] = operanda[5];
				operand[11][MEb] = operandb[5];
				operand[12][MEa] = operanda[6];
				operand[14][MEb] = operandb[6];
				operand[13][MEa] = operanda[7];
				operand[15][MEb] = operandb[7];
				operand[16][MEa] = operanda[8];
				operand[18][MEb] = operandb[8];
				operand[17][MEa] = operanda[9];
				operand[19][MEb] = operandb[9];
				operand[20][MEa] = operanda[10];
				operand[22][MEb] = operandb[10];
				operand[21][MEa] = operanda[11];
				operand[23][MEb] = operandb[11];
				operand[24][MEa] = operanda[12];
				operand[26][MEb] = operandb[12];
				operand[25][MEa] = operanda[13];
				operand[27][MEb] = operandb[13];
				operand[28][MEa] = operanda[14];
				operand[30][MEb] = operandb[14];
				operand[29][MEa] = operanda[15];
				operand[31][MEb] = operandb[15];
			}
			else if (i == 2)
			{
				operand[0][MEa] = operanda[0];
				operand[4][MEb] = operandb[0];
				operand[1][MEa] = operanda[1];
				operand[5][MEb] = operandb[1];
				operand[2][MEa] = operanda[2];
				operand[6][MEb] = operandb[2];
				operand[3][MEa] = operanda[3];
				operand[7][MEb] = operandb[3];
				operand[8][MEa] = operanda[4];
				operand[12][MEb] = operandb[4];
				operand[9][MEa] = operanda[5];
				operand[13][MEb] = operandb[5];
				operand[10][MEa] = operanda[6];
				operand[14][MEb] = operandb[6];
				operand[11][MEa] = operanda[7];
				operand[15][MEb] = operandb[7];
				operand[16][MEa] = operanda[8];
				operand[20][MEb] = operandb[8];
				operand[17][MEa] = operanda[9];
				operand[21][MEb] = operandb[9];
				operand[18][MEa] = operanda[10];
				operand[22][MEb] = operandb[10];
				operand[19][MEa] = operanda[11];
				operand[23][MEb] = operandb[11];
				operand[24][MEa] = operanda[12];
				operand[28][MEb] = operandb[12];
				operand[25][MEa] = operanda[13];
				operand[29][MEb] = operandb[13];
				operand[26][MEa] = operanda[14];
				operand[30][MEb] = operandb[14];
				operand[27][MEa] = operanda[15];
				operand[31][MEb] = operandb[15];
			}
			else if (i == 3)
			{
				operand[0][MEa] = operanda[0];
				operand[8][MEb] = operandb[0];
				operand[1][MEa] = operanda[1];
				operand[9][MEb] = operandb[1];
				operand[2][MEa] = operanda[2];
				operand[10][MEb] = operandb[2];
				operand[3][MEa] = operanda[3];
				operand[11][MEb] = operandb[3];
				operand[4][MEa] = operanda[4];
				operand[12][MEb] = operandb[4];
				operand[5][MEa] = operanda[5];
				operand[13][MEb] = operandb[5];
				operand[6][MEa] = operanda[6];
				operand[14][MEb] = operandb[6];
				operand[7][MEa] = operanda[7];
				operand[15][MEb] = operandb[7];
				operand[16][MEa] = operanda[8];
				operand[24][MEb] = operandb[8];
				operand[17][MEa] = operanda[9];
				operand[25][MEb] = operandb[9];
				operand[18][MEa] = operanda[10];
				operand[26][MEb] = operandb[10];
				operand[19][MEa] = operanda[11];
				operand[27][MEb] = operandb[11];
				operand[20][MEa] = operanda[12];
				operand[28][MEb] = operandb[12];
				operand[21][MEa] = operanda[13];
				operand[29][MEb] = operandb[13];
				operand[22][MEa] = operanda[14];
				operand[30][MEb] = operandb[14];
				operand[23][MEa] = operanda[15];
				operand[31][MEb] = operandb[15];
			}
			else if (i == 4)
			{
				operand[0][MEa] = operanda[0];
				operand[16][MEb] = operandb[0];
				operand[1][MEa] = operanda[1];
				operand[17][MEb] = operandb[1];
				operand[2][MEa] = operanda[2];
				operand[18][MEb] = operandb[2];
				operand[3][MEa] = operanda[3];
				operand[19][MEb] = operandb[3];
				operand[4][MEa] = operanda[4];
				operand[20][MEb] = operandb[4];
				operand[5][MEa] = operanda[5];
				operand[21][MEb] = operandb[5];
				operand[6][MEa] = operanda[6];
				operand[22][MEb] = operandb[6];
				operand[7][MEa] = operanda[7];
				operand[23][MEb] = operandb[7];
				operand[8][MEa] = operanda[8];
				operand[24][MEb] = operandb[8];
				operand[9][MEa] = operanda[9];
				operand[25][MEb] = operandb[9];
				operand[10][MEa] = operanda[10];
				operand[26][MEb] = operandb[10];
				operand[11][MEa] = operanda[11];
				operand[27][MEb] = operandb[11];
				operand[12][MEa] = operanda[12];
				operand[28][MEb] = operandb[12];
				operand[13][MEa] = operanda[13];
				operand[29][MEb] = operandb[13];
				operand[14][MEa] = operanda[14];
				operand[30][MEb] = operandb[14];
				operand[15][MEa] = operanda[15];
				operand[31][MEb] = operandb[15];
			}

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1:
	for (uint32_t i = stage_num - stage1_max; i < stage_num; i++)
	{
		uint32_t stepsize = step_num >> (i + 1);

		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			UDTYPE w0 = idt[temp1 + (j >> temp2)];
			UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operanda[m] = operand[m][MEa];
				operandb[m] = operand[m][MEb];
				operanda_[m] = operand[m + corenum][MEa];
				operandb_[m] = operand[m + corenum][MEb];
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m]);
			}
			j++;
			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m]);
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operand[m][MEa] = operanda[m];
				operand[m][MEb] = operandb[m];
				operand[m + corenum][MEa] = operanda_[m];
				operand[m + corenum][MEb] = operandb_[m];
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}
}

template <unsigned modcount, unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax>
void intt_1core_mods(UDTYPE operand[modcount][bramnum][bramsize],
					 UDTYPE idt[N], UDTYPE sidt[N], UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				uint32_t x = (arrayWindex[m] + temp_adder);
				//				w[m]=rp[x];
				//				sw[m]=srp[x];
				w[m] = idt[0];
				sw[m] = sidt[0];
			}

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				operanda[s][0] = operand[s][0][MEa];
				operandb[s][0] = operand[s][1][MEb];
			}

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w[0], sw[0], modulus, two_times_modulus, operanda[s][0], operandb[s][0]);
			}

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				operand[s][0][MEa] = operanda[s][0];
				operand[s][1][MEb] = operandb[s][0];
			}

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1:
	for (uint32_t i = stage_num - stage1_max; i < stage_num; i++)
	{
		uint32_t stepsize = step_num >> (i + 1);

		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			UDTYPE w0 = idt[temp1 + (j >> temp2)];
			UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				operanda[s][0] = operand[s][0][MEa];
				operandb[s][0] = operand[s][0][MEb];
				operanda_[s][0] = operand[s][0 + corenum][MEa];
				operandb_[s][0] = operand[s][0 + corenum][MEb];
			}

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[s][0], operandb[s][0]);
			}

			j++;

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				inverse_ntt_negacyclic_harvey_lazy_core(w00, sw00, modulus, two_times_modulus, operanda_[s][0], operandb_[s][0]);
			}

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				operand[s][0][MEa] = operanda[s][0];
				operand[s][0][MEb] = operandb[s][0];
				operand[s][0 + corenum][MEa] = operanda_[s][0];
				operand[s][0 + corenum][MEb] = operandb_[s][0];
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}
}

template <unsigned modcount, unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax>
void intt_2core_mods(UDTYPE operand[modcount][bramnum][bramsize],
					 UDTYPE idt[N], UDTYPE sidt[N], UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				uint32_t x = (arrayWindex[m] + temp_adder);
				//				w[m]=rp[x];
				//				sw[m]=srp[x];
				w[m] = idt[0];
				sw[m] = sidt[0];
			}

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				if (i == 0)
				{
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][1][MEb];
					operanda[s][1] = operand[s][2][MEa];
					operandb[s][1] = operand[s][3][MEb];
				}
				else if (i == 1)
				{
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][2][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][3][MEb];
				}
			}

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					inverse_ntt_negacyclic_harvey_lazy_core(w[m], sw[m], modulus, two_times_modulus, operanda[s][m], operandb[s][m]);
				}
			}

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				if (i == 0)
				{
					operand[s][0][MEa] = operanda[s][0];
					operand[s][1][MEb] = operandb[s][0];
					operand[s][2][MEa] = operanda[s][1];
					operand[s][3][MEb] = operandb[s][1];
				}
				else if (i == 1)
				{
					operand[s][0][MEa] = operanda[s][0];
					operand[s][2][MEb] = operandb[s][0];
					operand[s][1][MEa] = operanda[s][1];
					operand[s][3][MEb] = operandb[s][1];
				}
			}

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1:
	for (uint32_t i = stage_num - stage1_max; i < stage_num; i++)
	{
		uint32_t stepsize = step_num >> (i + 1);

		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			UDTYPE w0 = idt[temp1 + (j >> temp2)];
			UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					operanda[s][m] = operand[s][m][MEa];
					operandb[s][m] = operand[s][m][MEb];
					operanda_[s][m] = operand[s][m + corenum][MEa];
					operandb_[s][m] = operand[s][m + corenum][MEb];
				}
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					inverse_ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[s][m], operandb[s][m]);
				}
			}
			j++;
			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					inverse_ntt_negacyclic_harvey_lazy_core(w00, sw00, modulus, two_times_modulus, operanda_[s][m], operandb_[s][m]);
				}
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					operand[s][m][MEa] = operanda[s][m];
					operand[s][m][MEb] = operandb[s][m];
					operand[s][m + corenum][MEa] = operanda_[s][m];
					operand[s][m + corenum][MEb] = operandb_[s][m];
				}
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}
}

template <unsigned modcount, unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax>
void intt_4core_mods(UDTYPE operand[modcount][bramnum][bramsize],
					 UDTYPE idt[N], UDTYPE sidt[N], UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE operanda_out[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_out complete dim = 0
	UDTYPE operandb_out[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_out complete dim = 0
	UDTYPE operanda_2out[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_2out complete dim = 0
	UDTYPE operandb_2out[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_2out complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				uint32_t x = (arrayWindex[m] + temp_adder);
				//				w[m]=rp[x];
				//				sw[m]=srp[x];
				w[m] = idt[0];
				sw[m] = sidt[0];
			}

			if (i == 0)
			{
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][1][MEb];
					operanda[s][1] = operand[s][2][MEa];
					operandb[s][1] = operand[s][3][MEb];
					operanda[s][2] = operand[s][4][MEa];
					operandb[s][2] = operand[s][5][MEb];
					operanda[s][3] = operand[s][6][MEa];
					operandb[s][3] = operand[s][7][MEb];

					for (uint32_t m = 0; m < corenum; m++)
					{
#pragma HLS UNROLL
						inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus, operanda[s][m], operandb[s][m],
										 operanda_[s][m], operandb_[s][m]);
					}

					operand[s][0][MEa] = operanda_[s][0];
					operand[s][1][MEb] = operandb_[s][0];
					operand[s][2][MEa] = operanda_[s][1];
					operand[s][3][MEb] = operandb_[s][1];
					operand[s][4][MEa] = operanda_[s][2];
					operand[s][5][MEb] = operandb_[s][2];
					operand[s][6][MEa] = operanda_[s][3];
					operand[s][7][MEb] = operandb_[s][3];
				}
			}
			else if (i == 1)
			{
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][2][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][3][MEb];
					operanda[s][2] = operand[s][4][MEa];
					operandb[s][2] = operand[s][6][MEb];
					operanda[s][3] = operand[s][5][MEa];
					operandb[s][3] = operand[s][7][MEb];

					for (uint32_t m = 0; m < corenum; m++)
					{
#pragma HLS UNROLL
						inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus, operanda[s][m], operandb[s][m],
										 operanda_[s][m], operandb_[s][m]);
					}

					operand[s][0][MEa] = operanda_[s][0];
					operand[s][2][MEb] = operandb_[s][0];
					operand[s][1][MEa] = operanda_[s][1];
					operand[s][3][MEb] = operandb_[s][1];
					operand[s][4][MEa] = operanda_[s][2];
					operand[s][6][MEb] = operandb_[s][2];
					operand[s][5][MEa] = operanda_[s][3];
					operand[s][7][MEb] = operandb_[s][3];
				}
			}
			else if (i == 2)
			{
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][4][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][5][MEb];
					operanda[s][2] = operand[s][2][MEa];
					operandb[s][2] = operand[s][6][MEb];
					operanda[s][3] = operand[s][3][MEa];
					operandb[s][3] = operand[s][7][MEb];

					for (uint32_t m = 0; m < corenum; m++)
					{
#pragma HLS UNROLL
						inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus, operanda[s][m], operandb[s][m],
										 operanda_[s][m], operandb_[s][m]);
					}

					operand[s][0][MEa] = operanda_[s][0];
					operand[s][4][MEb] = operandb_[s][0];
					operand[s][1][MEa] = operanda_[s][1];
					operand[s][5][MEb] = operandb_[s][1];
					operand[s][2][MEa] = operanda_[s][2];
					operand[s][6][MEb] = operandb_[s][2];
					operand[s][3][MEa] = operanda_[s][3];
					operand[s][7][MEb] = operandb_[s][3];
				}
			}

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1:
	for (uint32_t i = stage_num - stage1_max; i < stage_num; i++)
	{
		// uint32_t stepsize = step_num >> (i + 1);

		uint32_t stepsize = 1 << (i - (stage_num - stage1_max));

		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			UDTYPE w0 = idt[temp1 + (j >> temp2)];
			UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					operanda[s][m] = operand[s][m][MEa];
					operandb[s][m] = operand[s][m][MEb];
					operanda_[s][m] = operand[s][m + corenum][MEa];
					operandb_[s][m] = operand[s][m + corenum][MEb];
				}
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[s][m], operandb[s][m],
									 operanda_out[s][m], operandb_out[s][m]);

					inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[s][m], operandb_[s][m],
									 operanda_2out[s][m], operandb_2out[s][m]);
				}
			}
			j++;
			/*for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[s][m], operandb_[s][m],
							operanda_2out[s][m], operandb_2out[s][m]);
				}
			}*/

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					/*std::cout << " m = " <<m << std::endl;
					std::cout << " m + corenum = " <<m + corenum << std::endl;
					std::cout << " MEa = " <<MEa << std::endl;
					std::cout << " MEb = " <<MEb << std::endl;*/
					operand[s][m][MEa] = operanda_out[s][m];
					operand[s][m][MEb] = operandb_out[s][m];
					operand[s][m + corenum][MEa] = operanda_2out[s][m];
					operand[s][m + corenum][MEb] = operandb_2out[s][m];
				}
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}
}

template <unsigned modcount, unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax>
void intt_8core_mods(UDTYPE operand[modcount][bramnum][bramsize],
					 UDTYPE idt[N], UDTYPE sidt[N], UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				uint32_t x = (arrayWindex[m] + temp_adder);
				//				w[m]=rp[x];
				//				sw[m]=srp[x];
				w[m] = idt[0];
				sw[m] = sidt[0];
			}

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				if (i == 0)
				{
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][1][MEb];
					operanda[s][1] = operand[s][2][MEa];
					operandb[s][1] = operand[s][3][MEb];
					operanda[s][2] = operand[s][4][MEa];
					operandb[s][2] = operand[s][5][MEb];
					operanda[s][3] = operand[s][6][MEa];
					operandb[s][3] = operand[s][7][MEb];
					operanda[s][4] = operand[s][8][MEa];
					operandb[s][4] = operand[s][9][MEb];
					operanda[s][5] = operand[s][10][MEa];
					operandb[s][5] = operand[s][11][MEb];
					operanda[s][6] = operand[s][12][MEa];
					operandb[s][6] = operand[s][13][MEb];
					operanda[s][7] = operand[s][14][MEa];
					operandb[s][7] = operand[s][15][MEb];
				}
				else if (i == 1)
				{
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][2][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][3][MEb];
					operanda[s][2] = operand[s][4][MEa];
					operandb[s][2] = operand[s][6][MEb];
					operanda[s][3] = operand[s][5][MEa];
					operandb[s][3] = operand[s][7][MEb];
					operanda[s][4] = operand[s][8][MEa];
					operandb[s][4] = operand[s][10][MEb];
					operanda[s][5] = operand[s][9][MEa];
					operandb[s][5] = operand[s][11][MEb];
					operanda[s][6] = operand[s][12][MEa];
					operandb[s][6] = operand[s][14][MEb];
					operanda[s][7] = operand[s][13][MEa];
					operandb[s][7] = operand[s][15][MEb];
				}
				else if (i == 2)
				{
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][4][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][5][MEb];
					operanda[s][2] = operand[s][2][MEa];
					operandb[s][2] = operand[s][6][MEb];
					operanda[s][3] = operand[s][3][MEa];
					operandb[s][3] = operand[s][7][MEb];
					operanda[s][4] = operand[s][8][MEa];
					operandb[s][4] = operand[s][12][MEb];
					operanda[s][5] = operand[s][9][MEa];
					operandb[s][5] = operand[s][13][MEb];
					operanda[s][6] = operand[s][10][MEa];
					operandb[s][6] = operand[s][14][MEb];
					operanda[s][7] = operand[s][11][MEa];
					operandb[s][7] = operand[s][15][MEb];
				}
				else if (i == 3)
				{
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][8][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][9][MEb];
					operanda[s][2] = operand[s][2][MEa];
					operandb[s][2] = operand[s][10][MEb];
					operanda[s][3] = operand[s][3][MEa];
					operandb[s][3] = operand[s][11][MEb];
					operanda[s][4] = operand[s][4][MEa];
					operandb[s][4] = operand[s][12][MEb];
					operanda[s][5] = operand[s][5][MEa];
					operandb[s][5] = operand[s][13][MEb];
					operanda[s][6] = operand[s][6][MEa];
					operandb[s][6] = operand[s][14][MEb];
					operanda[s][7] = operand[s][7][MEa];
					operandb[s][7] = operand[s][15][MEb];
				}
			}

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					inverse_ntt_negacyclic_harvey_lazy_core(w[m], sw[m], modulus, two_times_modulus, operanda[s][m], operandb[s][m]);
				}
			}

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				if (i == 0)
				{
					operand[s][0][MEa] = operanda[s][0];
					operand[s][1][MEb] = operandb[s][0];
					operand[s][2][MEa] = operanda[s][1];
					operand[s][3][MEb] = operandb[s][1];
					operand[s][4][MEa] = operanda[s][2];
					operand[s][5][MEb] = operandb[s][2];
					operand[s][6][MEa] = operanda[s][3];
					operand[s][7][MEb] = operandb[s][3];
					operand[s][8][MEa] = operanda[s][4];
					operand[s][9][MEb] = operandb[s][4];
					operand[s][10][MEa] = operanda[s][5];
					operand[s][11][MEb] = operandb[s][5];
					operand[s][12][MEa] = operanda[s][6];
					operand[s][13][MEb] = operandb[s][6];
					operand[s][14][MEa] = operanda[s][7];
					operand[s][15][MEb] = operandb[s][7];
				}
				else if (i == 1)
				{
					operand[s][0][MEa] = operanda[s][0];
					operand[s][2][MEb] = operandb[s][0];
					operand[s][1][MEa] = operanda[s][1];
					operand[s][3][MEb] = operandb[s][1];
					operand[s][4][MEa] = operanda[s][2];
					operand[s][6][MEb] = operandb[s][2];
					operand[s][5][MEa] = operanda[s][3];
					operand[s][7][MEb] = operandb[s][3];
					operand[s][8][MEa] = operanda[s][4];
					operand[s][10][MEb] = operandb[s][4];
					operand[s][9][MEa] = operanda[s][5];
					operand[s][11][MEb] = operandb[s][5];
					operand[s][12][MEa] = operanda[s][6];
					operand[s][14][MEb] = operandb[s][6];
					operand[s][13][MEa] = operanda[s][7];
					operand[s][15][MEb] = operandb[s][7];
				}
				else if (i == 2)
				{
					operand[s][0][MEa] = operanda[s][0];
					operand[s][4][MEb] = operandb[s][0];
					operand[s][1][MEa] = operanda[s][1];
					operand[s][5][MEb] = operandb[s][1];
					operand[s][2][MEa] = operanda[s][2];
					operand[s][6][MEb] = operandb[s][2];
					operand[s][3][MEa] = operanda[s][3];
					operand[s][7][MEb] = operandb[s][3];
					operand[s][8][MEa] = operanda[s][4];
					operand[s][12][MEb] = operandb[s][4];
					operand[s][9][MEa] = operanda[s][5];
					operand[s][13][MEb] = operandb[s][5];
					operand[s][10][MEa] = operanda[s][6];
					operand[s][14][MEb] = operandb[s][6];
					operand[s][11][MEa] = operanda[s][7];
					operand[s][15][MEb] = operandb[s][7];
				}
				else if (i == 3)
				{
					operand[s][0][MEa] = operanda[s][0];
					operand[s][8][MEb] = operandb[s][0];
					operand[s][1][MEa] = operanda[s][1];
					operand[s][9][MEb] = operandb[s][1];
					operand[s][2][MEa] = operanda[s][2];
					operand[s][10][MEb] = operandb[s][2];
					operand[s][3][MEa] = operanda[s][3];
					operand[s][11][MEb] = operandb[s][3];
					operand[s][4][MEa] = operanda[s][4];
					operand[s][12][MEb] = operandb[s][4];
					operand[s][5][MEa] = operanda[s][5];
					operand[s][13][MEb] = operandb[s][5];
					operand[s][6][MEa] = operanda[s][6];
					operand[s][14][MEb] = operandb[s][6];
					operand[s][7][MEa] = operanda[s][7];
					operand[s][15][MEb] = operandb[s][7];
				}
			}

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1:
	for (uint32_t i = stage_num - stage1_max; i < stage_num; i++)
	{
		uint32_t stepsize = step_num >> (i + 1);

		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			UDTYPE w0 = idt[temp1 + (j >> temp2)];
			UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					operanda[s][m] = operand[s][m][MEa];
					operandb[s][m] = operand[s][m][MEb];
					operanda_[s][m] = operand[s][m + corenum][MEa];
					operandb_[s][m] = operand[s][m + corenum][MEb];
				}
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					inverse_ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[s][m], operandb[s][m]);
				}
			}
			j++;
			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					inverse_ntt_negacyclic_harvey_lazy_core(w00, sw00, modulus, two_times_modulus, operanda_[s][m], operandb_[s][m]);
				}
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					operand[s][m][MEa] = operanda[s][m];
					operand[s][m][MEb] = operandb[s][m];
					operand[s][m + corenum][MEa] = operanda_[s][m];
					operand[s][m + corenum][MEb] = operandb_[s][m];
				}
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}
}

template <unsigned total_poly, unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax, unsigned type>
void intt_8core_mods_new(UDTYPE operand[total_poly][bramnum][bramsize],
						 UDTYPE idt[RPBRAMNUM][RPBRAMSIZE],
						 UDTYPE sidt[RPBRAMNUM][RPBRAMSIZE],
						 UDTYPE modulus)
{
#pragma HLS INLINE off

    UDTYPE operanda[total_poly][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
    UDTYPE operandb[total_poly][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
    UDTYPE operanda_[total_poly][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
    UDTYPE operandb_[total_poly][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

    UDTYPE operanda_out[total_poly][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_out complete dim = 0
    UDTYPE operandb_out[total_poly][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_out complete dim = 0
    UDTYPE operanda_2out[total_poly][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_2out complete dim = 0
    UDTYPE operandb_2out[total_poly][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_2out complete dim = 0

    UDTYPE two_times_modulus = modulus << 1;
    UDTYPE one_u64 = 1;
    uint32_t step_num = bramsize;
    uint32_t stage_num = STAGENUM;
    uint32_t stage1_max = stagemax;

    uint32_t MEa;
    uint32_t MEb;

    uint32_t arrayWindex;
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

    intt_stage2:
    for (uint32_t i = 0; i < stage_num - stage1_max; i++) {

        uint32_t temp4 = stage_num - 1 - i;
        uint32_t temp_par = corenum >> i;
        // uint32_t temp5 = i - stage1_max;

        MEa = 0;
        MEb = 0;

        arrayWindex = (1 << temp4);

        intt_stage2_inner:
        for (uint32_t j = 0; j < step_num; j++) {
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

            UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
            UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

            uint32_t temp_adder = temp_par * j;


            if (i == 0) {
                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    IDXTYPE x = (arrayWindex + temp_adder);
                    IDXTYPE rp_idx0_j;
                    rp_idx0_j = x >> logRPBRAMNUM;

                    w[0] = idt[0][rp_idx0_j];
                    sw[0] = sidt[0][rp_idx0_j];
                    w[1] = idt[1][rp_idx0_j];
                    sw[1] = sidt[1][rp_idx0_j];
                    w[2] = idt[2][rp_idx0_j];
                    sw[2] = sidt[2][rp_idx0_j];
                    w[3] = idt[3][rp_idx0_j];
                    sw[3] = sidt[3][rp_idx0_j];

                    w[4] = idt[0][rp_idx0_j];
                    w[5] = idt[1][rp_idx0_j];
                    w[6] = idt[2][rp_idx0_j];
                    w[7] = idt[3][rp_idx0_j];
                    sw[4] = sidt[0][rp_idx0_j];
                    sw[5] = sidt[1][rp_idx0_j];
                    sw[6] = sidt[2][rp_idx0_j];
                    sw[7] = sidt[3][rp_idx0_j];

                    operanda[s][0] = operand[s][0][MEa];
                    operandb[s][0] = operand[s][1][MEb];
                    operanda[s][1] = operand[s][2][MEa];
                    operandb[s][1] = operand[s][3][MEb];
                    operanda[s][2] = operand[s][4][MEa];
                    operandb[s][2] = operand[s][5][MEb];
                    operanda[s][3] = operand[s][6][MEa];
                    operandb[s][3] = operand[s][7][MEb];
                    operanda[s][4] = operand[s][8][MEa];
                    operandb[s][4] = operand[s][9][MEb];
                    operanda[s][5] = operand[s][10][MEa];
                    operandb[s][5] = operand[s][11][MEb];
                    operanda[s][6] = operand[s][12][MEa];
                    operandb[s][6] = operand[s][13][MEb];
                    operanda[s][7] = operand[s][14][MEa];
                    operandb[s][7] = operand[s][15][MEb];

                    for (uint32_t m = 0; m < corenum; m++) {
#pragma HLS UNROLL

                        inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
                                         operanda[s][m], operandb[s][m], operanda_[s][m], operandb_[s][m]);
                    }

                    operand[s][0][MEa] = operanda_[s][0];
                    operand[s][1][MEb] = operandb_[s][0];
                    operand[s][2][MEa] = operanda_[s][1];
                    operand[s][3][MEb] = operandb_[s][1];
                    operand[s][4][MEa] = operanda_[s][2];
                    operand[s][5][MEb] = operandb_[s][2];
                    operand[s][6][MEa] = operanda_[s][3];
                    operand[s][7][MEb] = operandb_[s][3];
                    operand[s][8][MEa] = operanda_[s][4];
                    operand[s][9][MEb] = operandb_[s][4];
                    operand[s][10][MEa] = operanda_[s][5];
                    operand[s][11][MEb] = operandb_[s][5];
                    operand[s][12][MEa] = operanda_[s][6];
                    operand[s][13][MEb] = operandb_[s][6];
                    operand[s][14][MEa] = operanda_[s][7];
                    operand[s][15][MEb] = operandb_[s][7];
                }
            }
            else if (i == 1) {
                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    IDXTYPE x = (arrayWindex + temp_adder);
                    IDXTYPE rp_idx0_i, rp_idx0_j;
                    rp_idx0_i = x(logRPBRAMNUM - 1, 0);
                    rp_idx0_j = x >> logRPBRAMNUM;

                    UDTYPE rp_t1, srp_t1;
                    UDTYPE rp_t2, srp_t2;

                    rp_t1 = idt[rp_idx0_i][rp_idx0_j];
                    srp_t1 = sidt[rp_idx0_i][rp_idx0_j];

                    rp_t2 = idt[rp_idx0_i + 1][rp_idx0_j];
                    srp_t2 = sidt[rp_idx0_i + 1][rp_idx0_j];

                    w[0] = rp_t1;
                    sw[0] = srp_t1;
                    w[1] = rp_t1;
                    sw[1] = srp_t1;
                    w[2] = rp_t1;
                    sw[2] = srp_t1;
                    w[3] = rp_t1;
                    sw[3] = srp_t1;

                    w[4] = rp_t2;
                    w[5] = rp_t2;
                    w[6] = rp_t2;
                    w[7] = rp_t2;
                    sw[4] = srp_t2;
                    sw[5] = srp_t2;
                    sw[6] = srp_t2;
                    sw[7] = srp_t2;

                    operanda[s][0] = operand[s][0][MEa];
                    operandb[s][0] = operand[s][2][MEb];
                    operanda[s][1] = operand[s][1][MEa];
                    operandb[s][1] = operand[s][3][MEb];
                    operanda[s][2] = operand[s][4][MEa];
                    operandb[s][2] = operand[s][6][MEb];
                    operanda[s][3] = operand[s][5][MEa];
                    operandb[s][3] = operand[s][7][MEb];
                    operanda[s][4] = operand[s][8][MEa];
                    operandb[s][4] = operand[s][10][MEb];
                    operanda[s][5] = operand[s][9][MEa];
                    operandb[s][5] = operand[s][11][MEb];
                    operanda[s][6] = operand[s][12][MEa];
                    operandb[s][6] = operand[s][14][MEb];
                    operanda[s][7] = operand[s][13][MEa];
                    operandb[s][7] = operand[s][15][MEb];

                    for (uint32_t m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
                                         operanda[s][m], operandb[s][m], operanda_[s][m], operandb_[s][m]);
                    }

                    operand[s][0][MEa] = operanda_[s][0];
                    operand[s][2][MEb] = operandb_[s][0];
                    operand[s][1][MEa] = operanda_[s][1];
                    operand[s][3][MEb] = operandb_[s][1];
                    operand[s][4][MEa] = operanda_[s][2];
                    operand[s][6][MEb] = operandb_[s][2];
                    operand[s][5][MEa] = operanda_[s][3];
                    operand[s][7][MEb] = operandb_[s][3];
                    operand[s][8][MEa] = operanda_[s][4];
                    operand[s][10][MEb] = operandb_[s][4];
                    operand[s][9][MEa] = operanda_[s][5];
                    operand[s][11][MEb] = operandb_[s][5];
                    operand[s][12][MEa] = operanda_[s][6];
                    operand[s][14][MEb] = operandb_[s][6];
                    operand[s][13][MEa] = operanda_[s][7];
                    operand[s][15][MEb] = operandb_[s][7];
                }
            }
            else if (i == 2) {
                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL

                    IDXTYPE x = (arrayWindex + temp_adder);
                    IDXTYPE rp_idx0_i, rp_idx0_j;
                    rp_idx0_i = x(logRPBRAMNUM - 1, 0);
                    rp_idx0_j = x >> logRPBRAMNUM;
                    // cout << "x = " << x << endl;

                    UDTYPE rp_t1, srp_t1;
                    UDTYPE rp_t2, srp_t2;
                    UDTYPE rp_t3, srp_t3;
                    UDTYPE rp_t4, srp_t4;

                    rp_t1 = idt[rp_idx0_i][rp_idx0_j];
                    srp_t1 = sidt[rp_idx0_i][rp_idx0_j];

                    rp_t2 = idt[rp_idx0_i + 1][rp_idx0_j];
                    srp_t2 = sidt[rp_idx0_i + 1][rp_idx0_j];

                    rp_t3 = idt[rp_idx0_i + 2][rp_idx0_j];
                    srp_t3 = sidt[rp_idx0_i + 2][rp_idx0_j];

                    rp_t4 = idt[rp_idx0_i + 3][rp_idx0_j];
                    srp_t4 = sidt[rp_idx0_i + 3][rp_idx0_j];

                    w[0] = rp_t1;
                    sw[0] = srp_t1;
                    w[1] = rp_t1;
                    sw[1] = srp_t1;

                    w[2] = rp_t2;
                    sw[2] = srp_t2;
                    w[3] = rp_t2;
                    sw[3] = srp_t2;

                    w[4] = rp_t3;
                    w[5] = rp_t3;
                    w[6] = rp_t4;
                    w[7] = rp_t4;
                    sw[4] = srp_t3;
                    sw[5] = srp_t3;
                    sw[6] = srp_t4;
                    sw[7] = srp_t4;

                    operanda[s][0] = operand[s][0][MEa];
                    operandb[s][0] = operand[s][4][MEb];
                    operanda[s][1] = operand[s][1][MEa];
                    operandb[s][1] = operand[s][5][MEb];
                    operanda[s][2] = operand[s][2][MEa];
                    operandb[s][2] = operand[s][6][MEb];
                    operanda[s][3] = operand[s][3][MEa];
                    operandb[s][3] = operand[s][7][MEb];
                    operanda[s][4] = operand[s][8][MEa];
                    operandb[s][4] = operand[s][12][MEb];
                    operanda[s][5] = operand[s][9][MEa];
                    operandb[s][5] = operand[s][13][MEb];
                    operanda[s][6] = operand[s][10][MEa];
                    operandb[s][6] = operand[s][14][MEb];
                    operanda[s][7] = operand[s][11][MEa];
                    operandb[s][7] = operand[s][15][MEb];

                    for (uint32_t m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
                                         operanda[s][m], operandb[s][m], operanda_[s][m], operandb_[s][m]);
                    }

                    operand[s][0][MEa] = operanda_[s][0];
                    operand[s][4][MEb] = operandb_[s][0];
                    operand[s][1][MEa] = operanda_[s][1];
                    operand[s][5][MEb] = operandb_[s][1];
                    operand[s][2][MEa] = operanda_[s][2];
                    operand[s][6][MEb] = operandb_[s][2];
                    operand[s][3][MEa] = operanda_[s][3];
                    operand[s][7][MEb] = operandb_[s][3];
                    operand[s][8][MEa] = operanda_[s][4];
                    operand[s][12][MEb] = operandb_[s][4];
                    operand[s][9][MEa] = operanda_[s][5];
                    operand[s][13][MEb] = operandb_[s][5];
                    operand[s][10][MEa] = operanda_[s][6];
                    operand[s][14][MEb] = operandb_[s][6];
                    operand[s][11][MEa] = operanda_[s][7];
                    operand[s][15][MEb] = operandb_[s][7];
                }
            }
            else if (i == 3) {
                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL

                    for (uint32_t m = 0; m < corenum; m++) {
#pragma HLS UNROLL

                        IDXTYPE x = (arrayWindex + temp_adder);
                        IDXTYPE rp_idx0_i, rp_idx0_j;
                        rp_idx0_i = x(logRPBRAMNUM - 1, 0);
                        rp_idx0_j = x >> logRPBRAMNUM;
                        // cout << "x = " << x << endl;

                        w[m] = idt[rp_idx0_i][rp_idx0_j];
                        sw[m] = sidt[rp_idx0_i][rp_idx0_j];
                    }

                    operanda[s][0] = operand[s][0][MEa];
                    operandb[s][0] = operand[s][8][MEb];
                    operanda[s][1] = operand[s][1][MEa];
                    operandb[s][1] = operand[s][9][MEb];
                    operanda[s][2] = operand[s][2][MEa];
                    operandb[s][2] = operand[s][10][MEb];
                    operanda[s][3] = operand[s][3][MEa];
                    operandb[s][3] = operand[s][11][MEb];
                    operanda[s][4] = operand[s][4][MEa];
                    operandb[s][4] = operand[s][12][MEb];
                    operanda[s][5] = operand[s][5][MEa];
                    operandb[s][5] = operand[s][13][MEb];
                    operanda[s][6] = operand[s][6][MEa];
                    operandb[s][6] = operand[s][14][MEb];
                    operanda[s][7] = operand[s][7][MEa];
                    operandb[s][7] = operand[s][15][MEb];

                    for (uint32_t m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        inverse_ntt_core(w[m], sw[m], modulus, two_times_modulus,
                                         operanda[s][m], operandb[s][m], operanda_[s][m], operandb_[s][m]);
                    }

                    operand[s][0][MEa] = operanda_[s][0];
                    operand[s][8][MEb] = operandb_[s][0];
                    operand[s][1][MEa] = operanda_[s][1];
                    operand[s][9][MEb] = operandb_[s][1];
                    operand[s][2][MEa] = operanda_[s][2];
                    operand[s][10][MEb] = operandb_[s][2];
                    operand[s][3][MEa] = operanda_[s][3];
                    operand[s][11][MEb] = operandb_[s][3];
                    operand[s][4][MEa] = operanda_[s][4];
                    operand[s][12][MEb] = operandb_[s][4];
                    operand[s][5][MEa] = operanda_[s][5];
                    operand[s][13][MEb] = operandb_[s][5];
                    operand[s][6][MEa] = operanda_[s][6];
                    operand[s][14][MEb] = operandb_[s][6];
                    operand[s][7][MEa] = operanda_[s][7];
                    operand[s][15][MEb] = operandb_[s][7];
                }
            }

            MEa += 1;
            MEb += 1;
        }
    }

    intt_stage1_without_last:
    for (uint32_t i = stage_num - stage1_max; i < stage_num - 1; i++) {

        uint32_t stepsize = 1 << (i - (stage_num - stage1_max));
        uint32_t temp2 = i - (stage_num - stage1_max - 1);
        uint32_t temp1 = step_num >> temp2;
        uint32_t temp3 = (one_u64 << temp2) - one_u64;
        uint32_t stage_tempa = 12;

        MEa = 0;
        MEb = stepsize;

        // cout << "stepsize" << stepsize << endl;

        intt_stage1_inner:
        for (uint32_t j = 0; j < step_num; j++) {
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

            // UDTYPE w0 = idt[temp1 + (j >> temp2)];
            // UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
            // UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
            // UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

            IDXTYPE rp_idx0 = temp1 + (j >> temp2);
            IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
            IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

            rp_idx0_i = rp_idx0(logRPBRAMNUM - 1, 0);
            rp_idx0_j = rp_idx0 >> logRPBRAMNUM;
            rp_idx00_i = rp_idx00(logRPBRAMNUM - 1, 0);
            rp_idx00_j = rp_idx00 >> logRPBRAMNUM;

            UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
            UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
            UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
            UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

            for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                    operanda[s][m] = operand[s][m][MEa];
                    operandb[s][m] = operand[s][m][MEb];
                    operanda_[s][m] = operand[s][m + corenum][MEa];
                    operandb_[s][m] = operand[s][m + corenum][MEb];
                }
            }

            for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL

                    inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[s][m], operandb[s][m],
                                     operanda_out[s][m], operandb_out[s][m]);

                    inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[s][m], operandb_[s][m],
                                     operanda_2out[s][m], operandb_2out[s][m]);
                }
            }
            j++;

            for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL

                    operand[s][m][MEa] = operanda_out[s][m];
                    operand[s][m][MEb] = operandb_out[s][m];
                    operand[s][m + corenum][MEa] = operanda_2out[s][m];
                    operand[s][m + corenum][MEb] = operandb_2out[s][m];
                }
            }

            if (((j + 1) & temp3) == 0) {
                MEa += stepsize + 1;
                MEb += stepsize + 1;
            }
            else {
                MEa += 1;
                MEb += 1;
            }
        }
    }

    if (type == 0) {
        type0_intt_stage1_last:
        for (uint32_t i = stage_num - 1; i < stage_num; i++) {

            uint32_t stepsize = 1 << (i - (stage_num - stage1_max));
            uint32_t temp2 = i - (stage_num - stage1_max - 1);
            uint32_t temp1 = step_num >> temp2;
            uint32_t temp3 = (one_u64 << temp2) - one_u64;
            uint32_t stage_tempa = 12;

            MEa = 0;
            MEb = stepsize;

            // cout << "stepsize" << stepsize << endl;

            type0_intt_stage1_inner_last:
            for (uint32_t j = 0; j < step_num; j++) {
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

                // UDTYPE w0 = idt[temp1 + (j >> temp2)];
                // UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
                // UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
                // UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

                IDXTYPE rp_idx0 = temp1 + (j >> temp2);
                IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
                IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

                rp_idx0_i = rp_idx0(logRPBRAMNUM - 1, 0);
                rp_idx0_j = rp_idx0 >> logRPBRAMNUM;
                rp_idx00_i = rp_idx00(logRPBRAMNUM - 1, 0);
                rp_idx00_j = rp_idx00 >> logRPBRAMNUM;

                UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
                UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
                UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
                UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        operanda[s][m] = operand[s][m][MEa];
                        operandb[s][m] = operand[s][m][MEb];
                        operanda_[s][m] = operand[s][m + corenum][MEa];
                        operandb_[s][m] = operand[s][m + corenum][MEb];
                    }
                }

                for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                    for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                        inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[s][m], operandb[s][m],
                                         operanda_out[s][m], operandb_out[s][m]);

                        inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[s][m], operandb_[s][m],
                                         operanda_2out[s][m], operandb_2out[s][m]);
                    }
                }
                j++;

                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        operand[s][m][MEa] = operanda_out[s][m];
                        operand[s][m][MEb] = operandb_out[s][m];
                        operand[s][m + corenum][MEa] = operanda_2out[s][m];
                        operand[s][m + corenum][MEb] = operandb_2out[s][m];
                    }
                }

                if (((j + 1) & temp3) == 0) {
                    MEa += stepsize + 1;
                    MEb += stepsize + 1;
                }
                else {
                    MEa += 1;
                    MEb += 1;
                }
            }
        }
    }
    else if (type == 1) {
        type1_intt_stage1_last:
        for (uint32_t i = stage_num - 1; i < stage_num; i++) {

            uint32_t stepsize = 1 << (i - (stage_num - stage1_max));
            uint32_t temp2 = i - (stage_num - stage1_max - 1);
            uint32_t temp1 = step_num >> temp2;
            uint32_t temp3 = (one_u64 << temp2) - one_u64;
            uint32_t stage_tempa = 12;

            MEa = 0;
            MEb = stepsize;

            // cout << "stepsize" << stepsize << endl;

            type1_intt_stage1_inner_last:
            for (uint32_t j = 0; j < step_num; j++) {
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

                IDXTYPE rp_idx0 = temp1 + (j >> temp2);
                IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
                IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

                rp_idx0_i = rp_idx0(logRPBRAMNUM - 1, 0);
                rp_idx0_j = rp_idx0 >> logRPBRAMNUM;
                rp_idx00_i = rp_idx00(logRPBRAMNUM - 1, 0);
                rp_idx00_j = rp_idx00 >> logRPBRAMNUM;

                UDTYPE w0 = idt[rp_idx0_i][rp_idx0_j];
                UDTYPE sw0 = sidt[rp_idx0_i][rp_idx0_j];
                UDTYPE w00 = idt[rp_idx00_i][rp_idx00_j];
                UDTYPE sw00 = sidt[rp_idx00_i][rp_idx00_j];

                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        operanda[s][m] = operand[s][m][MEa];
                        operandb[s][m] = operand[s][m][MEb];
                        operanda_[s][m] = operand[s][m + corenum][MEa];
                        operandb_[s][m] = operand[s][m + corenum][MEb];
                    }
                }

                for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                    for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                        inverse_ntt_core(w0, sw0, modulus, two_times_modulus, operanda[s][m], operandb[s][m],
                                         operanda_out[s][m], operandb_out[s][m]);

                        inverse_ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[s][m], operandb_[s][m],
                                         operanda_2out[s][m], operandb_2out[s][m]);
                    }
                }
                j++;

                for (int s = 0; s < total_poly; s++) {
#pragma HLS UNROLL
                    for (int m = 0; m < corenum; m++) {
#pragma HLS UNROLL
                        operand[s][m][MEa] = static_cast<UDTYPE>(operanda_out[s][m]) - (modulus &
                                                                                        static_cast<UDTYPE>(-static_cast<DTYPE>(
                                                                                                operanda_out[s][m] >=
                                                                                                modulus)));
                        operand[s][m][MEb] = static_cast<UDTYPE>(operandb_out[s][m]) - (modulus &
                                                                                        static_cast<UDTYPE>(-static_cast<DTYPE>(
                                                                                                operandb_out[s][m] >=
                                                                                                modulus)));
                        operand[s][m + corenum][MEa] = static_cast<UDTYPE>(operanda_2out[s][m]) - (modulus &
                                                                                                   static_cast<UDTYPE>(-static_cast<DTYPE>(
                                                                                                           operanda_2out[s][m] >=
                                                                                                           modulus)));
                        operand[s][m + corenum][MEb] = static_cast<UDTYPE>(operandb_2out[s][m]) - (modulus &
                                                                                                   static_cast<UDTYPE>(-static_cast<DTYPE>(
                                                                                                           operandb_2out[s][m] >=
                                                                                                           modulus)));
                    }
                }

                if (((j + 1) & temp3) == 0) {
                    MEa += stepsize + 1;
                    MEb += stepsize + 1;
                }
                else {
                    MEa += 1;
                    MEb += 1;
                }
            }
        }
    }
}


template <unsigned modcount, unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax>
void intt_16core_mods(UDTYPE operand[modcount][bramnum][bramsize],
					  UDTYPE idt[N], UDTYPE sidt[N], UDTYPE modulus)
{
#pragma HLS INLINE off

	UDTYPE operanda[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[modcount][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

intt_stage2:
	for (uint32_t i = 0; i < stage_num - stage1_max; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> i;
		// uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << temp4) + (m >> i);
		}

	intt_stage2_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				uint32_t x = (arrayWindex[m] + temp_adder);
				//				w[m]=rp[x];
				//				sw[m]=srp[x];
				w[m] = idt[0];
				sw[m] = sidt[0];
			}

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				if (i == 0)
				{
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][1][MEb];
					operanda[s][1] = operand[s][2][MEa];
					operandb[s][1] = operand[s][3][MEb];
					operanda[s][2] = operand[s][4][MEa];
					operandb[s][2] = operand[s][5][MEb];
					operanda[s][3] = operand[s][6][MEa];
					operandb[s][3] = operand[s][7][MEb];
					operanda[s][4] = operand[s][8][MEa];
					operandb[s][4] = operand[s][9][MEb];
					operanda[s][5] = operand[s][10][MEa];
					operandb[s][5] = operand[s][11][MEb];
					operanda[s][6] = operand[s][12][MEa];
					operandb[s][6] = operand[s][13][MEb];
					operanda[s][7] = operand[s][14][MEa];
					operandb[s][7] = operand[s][15][MEb];
					operanda[s][8] = operand[s][16][MEa];
					operandb[s][8] = operand[s][17][MEb];
					operanda[s][9] = operand[s][18][MEa];
					operandb[s][9] = operand[s][19][MEb];
					operanda[s][10] = operand[s][20][MEa];
					operandb[s][10] = operand[s][21][MEb];
					operanda[s][11] = operand[s][22][MEa];
					operandb[s][11] = operand[s][23][MEb];
					operanda[s][12] = operand[s][24][MEa];
					operandb[s][12] = operand[s][25][MEb];
					operanda[s][13] = operand[s][26][MEa];
					operandb[s][13] = operand[s][27][MEb];
					operanda[s][14] = operand[s][28][MEa];
					operandb[s][14] = operand[s][29][MEb];
					operanda[s][15] = operand[s][30][MEa];
					operandb[s][15] = operand[s][31][MEb];
				}
				else if (i == 1)
				{
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][2][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][3][MEb];
					operanda[s][2] = operand[s][4][MEa];
					operandb[s][2] = operand[s][6][MEb];
					operanda[s][3] = operand[s][5][MEa];
					operandb[s][3] = operand[s][7][MEb];
					operanda[s][4] = operand[s][8][MEa];
					operandb[s][4] = operand[s][10][MEb];
					operanda[s][5] = operand[s][9][MEa];
					operandb[s][5] = operand[s][11][MEb];
					operanda[s][6] = operand[s][12][MEa];
					operandb[s][6] = operand[s][14][MEb];
					operanda[s][7] = operand[s][13][MEa];
					operandb[s][7] = operand[s][15][MEb];
					operanda[s][8] = operand[s][16][MEa];
					operandb[s][8] = operand[s][18][MEb];
					operanda[s][9] = operand[s][17][MEa];
					operandb[s][9] = operand[s][19][MEb];
					operanda[s][10] = operand[s][20][MEa];
					operandb[s][10] = operand[s][22][MEb];
					operanda[s][11] = operand[s][21][MEa];
					operandb[s][11] = operand[s][23][MEb];
					operanda[s][12] = operand[s][24][MEa];
					operandb[s][12] = operand[s][26][MEb];
					operanda[s][13] = operand[s][25][MEa];
					operandb[s][13] = operand[s][27][MEb];
					operanda[s][14] = operand[s][28][MEa];
					operandb[s][14] = operand[s][30][MEb];
					operanda[s][15] = operand[s][29][MEa];
					operandb[s][15] = operand[s][31][MEb];
				}
				else if (i == 2)
				{
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][4][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][5][MEb];
					operanda[s][2] = operand[s][2][MEa];
					operandb[s][2] = operand[s][6][MEb];
					operanda[s][3] = operand[s][3][MEa];
					operandb[s][3] = operand[s][7][MEb];
					operanda[s][4] = operand[s][8][MEa];
					operandb[s][4] = operand[s][12][MEb];
					operanda[s][5] = operand[s][9][MEa];
					operandb[s][5] = operand[s][13][MEb];
					operanda[s][6] = operand[s][10][MEa];
					operandb[s][6] = operand[s][14][MEb];
					operanda[s][7] = operand[s][11][MEa];
					operandb[s][7] = operand[s][15][MEb];
					operanda[s][8] = operand[s][16][MEa];
					operandb[s][8] = operand[s][20][MEb];
					operanda[s][9] = operand[s][17][MEa];
					operandb[s][9] = operand[s][21][MEb];
					operanda[s][10] = operand[s][18][MEa];
					operandb[s][10] = operand[s][22][MEb];
					operanda[s][11] = operand[s][19][MEa];
					operandb[s][11] = operand[s][23][MEb];
					operanda[s][12] = operand[s][24][MEa];
					operandb[s][12] = operand[s][28][MEb];
					operanda[s][13] = operand[s][25][MEa];
					operandb[s][13] = operand[s][29][MEb];
					operanda[s][14] = operand[s][26][MEa];
					operandb[s][14] = operand[s][30][MEb];
					operanda[s][15] = operand[s][27][MEa];
					operandb[s][15] = operand[s][31][MEb];
				}
				else if (i == 3)
				{
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][8][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][9][MEb];
					operanda[s][2] = operand[s][2][MEa];
					operandb[s][2] = operand[s][10][MEb];
					operanda[s][3] = operand[s][3][MEa];
					operandb[s][3] = operand[s][11][MEb];
					operanda[s][4] = operand[s][4][MEa];
					operandb[s][4] = operand[s][12][MEb];
					operanda[s][5] = operand[s][5][MEa];
					operandb[s][5] = operand[s][13][MEb];
					operanda[s][6] = operand[s][6][MEa];
					operandb[s][6] = operand[s][14][MEb];
					operanda[s][7] = operand[s][7][MEa];
					operandb[s][7] = operand[s][15][MEb];
					operanda[s][8] = operand[s][16][MEa];
					operandb[s][8] = operand[s][24][MEb];
					operanda[s][9] = operand[s][17][MEa];
					operandb[s][9] = operand[s][25][MEb];
					operanda[s][10] = operand[s][18][MEa];
					operandb[s][10] = operand[s][26][MEb];
					operanda[s][11] = operand[s][19][MEa];
					operandb[s][11] = operand[s][27][MEb];
					operanda[s][12] = operand[s][20][MEa];
					operandb[s][12] = operand[s][28][MEb];
					operanda[s][13] = operand[s][21][MEa];
					operandb[s][13] = operand[s][29][MEb];
					operanda[s][14] = operand[s][22][MEa];
					operandb[s][14] = operand[s][30][MEb];
					operanda[s][15] = operand[s][23][MEa];
					operandb[s][15] = operand[s][31][MEb];
				}
				else if (i == 4)
				{
					operanda[s][0] = operand[s][0][MEa];
					operandb[s][0] = operand[s][16][MEb];
					operanda[s][1] = operand[s][1][MEa];
					operandb[s][1] = operand[s][17][MEb];
					operanda[s][2] = operand[s][2][MEa];
					operandb[s][2] = operand[s][18][MEb];
					operanda[s][3] = operand[s][3][MEa];
					operandb[s][3] = operand[s][19][MEb];
					operanda[s][4] = operand[s][4][MEa];
					operandb[s][4] = operand[s][20][MEb];
					operanda[s][5] = operand[s][5][MEa];
					operandb[s][5] = operand[s][21][MEb];
					operanda[s][6] = operand[s][6][MEa];
					operandb[s][6] = operand[s][22][MEb];
					operanda[s][7] = operand[s][7][MEa];
					operandb[s][7] = operand[s][23][MEb];
					operanda[s][8] = operand[s][8][MEa];
					operandb[s][8] = operand[s][24][MEb];
					operanda[s][9] = operand[s][9][MEa];
					operandb[s][9] = operand[s][25][MEb];
					operanda[s][10] = operand[s][10][MEa];
					operandb[s][10] = operand[s][26][MEb];
					operanda[s][11] = operand[s][11][MEa];
					operandb[s][11] = operand[s][27][MEb];
					operanda[s][12] = operand[s][12][MEa];
					operandb[s][12] = operand[s][28][MEb];
					operanda[s][13] = operand[s][13][MEa];
					operandb[s][13] = operand[s][29][MEb];
					operanda[s][14] = operand[s][14][MEa];
					operandb[s][14] = operand[s][30][MEb];
					operanda[s][15] = operand[s][15][MEa];
					operandb[s][15] = operand[s][31][MEb];
				}
			}

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					inverse_ntt_negacyclic_harvey_lazy_core(w[m], sw[m], modulus, two_times_modulus, operanda[s][m], operandb[s][m]);
				}
			}

			for (int s = 0; s < modcount; s++)
			{
#pragma HLS UNROLL
				if (i == 0)
				{
					operand[s][0][MEa] = operanda[s][0];
					operand[s][1][MEb] = operandb[s][0];
					operand[s][2][MEa] = operanda[s][1];
					operand[s][3][MEb] = operandb[s][1];
					operand[s][4][MEa] = operanda[s][2];
					operand[s][5][MEb] = operandb[s][2];
					operand[s][6][MEa] = operanda[s][3];
					operand[s][7][MEb] = operandb[s][3];
					operand[s][8][MEa] = operanda[s][4];
					operand[s][9][MEb] = operandb[s][4];
					operand[s][10][MEa] = operanda[s][5];
					operand[s][11][MEb] = operandb[s][5];
					operand[s][12][MEa] = operanda[s][6];
					operand[s][13][MEb] = operandb[s][6];
					operand[s][14][MEa] = operanda[s][7];
					operand[s][15][MEb] = operandb[s][7];
					operand[s][16][MEa] = operanda[s][8];
					operand[s][17][MEb] = operandb[s][8];
					operand[s][18][MEa] = operanda[s][9];
					operand[s][19][MEb] = operandb[s][9];
					operand[s][20][MEa] = operanda[s][10];
					operand[s][21][MEb] = operandb[s][10];
					operand[s][22][MEa] = operanda[s][11];
					operand[s][23][MEb] = operandb[s][11];
					operand[s][24][MEa] = operanda[s][12];
					operand[s][25][MEb] = operandb[s][12];
					operand[s][26][MEa] = operanda[s][13];
					operand[s][27][MEb] = operandb[s][13];
					operand[s][28][MEa] = operanda[s][14];
					operand[s][29][MEb] = operandb[s][14];
					operand[s][30][MEa] = operanda[s][15];
					operand[s][31][MEb] = operandb[s][15];
				}
				else if (i == 1)
				{
					operand[s][0][MEa] = operanda[s][0];
					operand[s][2][MEb] = operandb[s][0];
					operand[s][1][MEa] = operanda[s][1];
					operand[s][3][MEb] = operandb[s][1];
					operand[s][4][MEa] = operanda[s][2];
					operand[s][6][MEb] = operandb[s][2];
					operand[s][5][MEa] = operanda[s][3];
					operand[s][7][MEb] = operandb[s][3];
					operand[s][8][MEa] = operanda[s][4];
					operand[s][10][MEb] = operandb[s][4];
					operand[s][9][MEa] = operanda[s][5];
					operand[s][11][MEb] = operandb[s][5];
					operand[s][12][MEa] = operanda[s][6];
					operand[s][14][MEb] = operandb[s][6];
					operand[s][13][MEa] = operanda[s][7];
					operand[s][15][MEb] = operandb[s][7];
					operand[s][16][MEa] = operanda[s][8];
					operand[s][18][MEb] = operandb[s][8];
					operand[s][17][MEa] = operanda[s][9];
					operand[s][19][MEb] = operandb[s][9];
					operand[s][20][MEa] = operanda[s][10];
					operand[s][22][MEb] = operandb[s][10];
					operand[s][21][MEa] = operanda[s][11];
					operand[s][23][MEb] = operandb[s][11];
					operand[s][24][MEa] = operanda[s][12];
					operand[s][26][MEb] = operandb[s][12];
					operand[s][25][MEa] = operanda[s][13];
					operand[s][27][MEb] = operandb[s][13];
					operand[s][28][MEa] = operanda[s][14];
					operand[s][30][MEb] = operandb[s][14];
					operand[s][29][MEa] = operanda[s][15];
					operand[s][31][MEb] = operandb[s][15];
				}
				else if (i == 2)
				{
					operand[s][0][MEa] = operanda[s][0];
					operand[s][4][MEb] = operandb[s][0];
					operand[s][1][MEa] = operanda[s][1];
					operand[s][5][MEb] = operandb[s][1];
					operand[s][2][MEa] = operanda[s][2];
					operand[s][6][MEb] = operandb[s][2];
					operand[s][3][MEa] = operanda[s][3];
					operand[s][7][MEb] = operandb[s][3];
					operand[s][8][MEa] = operanda[s][4];
					operand[s][12][MEb] = operandb[s][4];
					operand[s][9][MEa] = operanda[s][5];
					operand[s][13][MEb] = operandb[s][5];
					operand[s][10][MEa] = operanda[s][6];
					operand[s][14][MEb] = operandb[s][6];
					operand[s][11][MEa] = operanda[s][7];
					operand[s][15][MEb] = operandb[s][7];
					operand[s][16][MEa] = operanda[s][8];
					operand[s][20][MEb] = operandb[s][8];
					operand[s][17][MEa] = operanda[s][9];
					operand[s][21][MEb] = operandb[s][9];
					operand[s][18][MEa] = operanda[s][10];
					operand[s][22][MEb] = operandb[s][10];
					operand[s][19][MEa] = operanda[s][11];
					operand[s][23][MEb] = operandb[s][11];
					operand[s][24][MEa] = operanda[s][12];
					operand[s][28][MEb] = operandb[s][12];
					operand[s][25][MEa] = operanda[s][13];
					operand[s][29][MEb] = operandb[s][13];
					operand[s][26][MEa] = operanda[s][14];
					operand[s][30][MEb] = operandb[s][14];
					operand[s][27][MEa] = operanda[s][15];
					operand[s][31][MEb] = operandb[s][15];
				}
				else if (i == 3)
				{
					operand[s][0][MEa] = operanda[s][0];
					operand[s][8][MEb] = operandb[s][0];
					operand[s][1][MEa] = operanda[s][1];
					operand[s][9][MEb] = operandb[s][1];
					operand[s][2][MEa] = operanda[s][2];
					operand[s][10][MEb] = operandb[s][2];
					operand[s][3][MEa] = operanda[s][3];
					operand[s][11][MEb] = operandb[s][3];
					operand[s][4][MEa] = operanda[s][4];
					operand[s][12][MEb] = operandb[s][4];
					operand[s][5][MEa] = operanda[s][5];
					operand[s][13][MEb] = operandb[s][5];
					operand[s][6][MEa] = operanda[s][6];
					operand[s][14][MEb] = operandb[s][6];
					operand[s][7][MEa] = operanda[s][7];
					operand[s][15][MEb] = operandb[s][7];
					operand[s][16][MEa] = operanda[s][8];
					operand[s][24][MEb] = operandb[s][8];
					operand[s][17][MEa] = operanda[s][9];
					operand[s][25][MEb] = operandb[s][9];
					operand[s][18][MEa] = operanda[s][10];
					operand[s][26][MEb] = operandb[s][10];
					operand[s][19][MEa] = operanda[s][11];
					operand[s][27][MEb] = operandb[s][11];
					operand[s][20][MEa] = operanda[s][12];
					operand[s][28][MEb] = operandb[s][12];
					operand[s][21][MEa] = operanda[s][13];
					operand[s][29][MEb] = operandb[s][13];
					operand[s][22][MEa] = operanda[s][14];
					operand[s][30][MEb] = operandb[s][14];
					operand[s][23][MEa] = operanda[s][15];
					operand[s][31][MEb] = operandb[s][15];
				}
				else if (i == 4)
				{
					operand[s][0][MEa] = operanda[s][0];
					operand[s][16][MEb] = operandb[s][0];
					operand[s][1][MEa] = operanda[s][1];
					operand[s][17][MEb] = operandb[s][1];
					operand[s][2][MEa] = operanda[s][2];
					operand[s][18][MEb] = operandb[s][2];
					operand[s][3][MEa] = operanda[s][3];
					operand[s][19][MEb] = operandb[s][3];
					operand[s][4][MEa] = operanda[s][4];
					operand[s][20][MEb] = operandb[s][4];
					operand[s][5][MEa] = operanda[s][5];
					operand[s][21][MEb] = operandb[s][5];
					operand[s][6][MEa] = operanda[s][6];
					operand[s][22][MEb] = operandb[s][6];
					operand[s][7][MEa] = operanda[s][7];
					operand[s][23][MEb] = operandb[s][7];
					operand[s][8][MEa] = operanda[s][8];
					operand[s][24][MEb] = operandb[s][8];
					operand[s][9][MEa] = operanda[s][9];
					operand[s][25][MEb] = operandb[s][9];
					operand[s][10][MEa] = operanda[s][10];
					operand[s][26][MEb] = operandb[s][10];
					operand[s][11][MEa] = operanda[s][11];
					operand[s][27][MEb] = operandb[s][11];
					operand[s][12][MEa] = operanda[s][12];
					operand[s][28][MEb] = operandb[s][12];
					operand[s][13][MEa] = operanda[s][13];
					operand[s][29][MEb] = operandb[s][13];
					operand[s][14][MEa] = operanda[s][14];
					operand[s][30][MEb] = operandb[s][14];
					operand[s][15][MEa] = operanda[s][15];
					operand[s][31][MEb] = operandb[s][15];
				}
			}

			MEa += 1;
			MEb += 1;
		}
	}

intt_stage1:
	for (uint32_t i = stage_num - stage1_max; i < stage_num; i++)
	{
		uint32_t stepsize = step_num >> (i + 1);

		uint32_t temp2 = i - (stage_num - stage1_max - 1);
		uint32_t temp1 = step_num >> temp2;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	intt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			UDTYPE w0 = idt[temp1 + (j >> temp2)];
			UDTYPE sw0 = sidt[temp1 + (j >> temp2)];
			UDTYPE w00 = idt[temp1 + ((j + 1) >> temp2)];
			UDTYPE sw00 = sidt[temp1 + ((j + 1) >> temp2)];

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					operanda[s][m] = operand[s][m][MEa];
					operandb[s][m] = operand[s][m][MEb];
					operanda_[s][m] = operand[s][m + corenum][MEa];
					operandb_[s][m] = operand[s][m + corenum][MEb];
				}
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					inverse_ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[s][m], operandb[s][m]);
				}
			}
			j++;

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					inverse_ntt_negacyclic_harvey_lazy_core(w00, sw00, modulus, two_times_modulus, operanda_[s][m], operandb_[s][m]);
				}
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				for (int s = 0; s < modcount; s++)
				{
#pragma HLS UNROLL
					operand[s][m][MEa] = operanda[s][m];
					operand[s][m][MEb] = operandb[s][m];
					operand[s][m + corenum][MEa] = operanda_[s][m];
					operand[s][m + corenum][MEb] = operandb_[s][m];
				}
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}
}

template <unsigned corenum, unsigned bramnum, unsigned bramsize>
void intt_double_1(UDTYPE operand1[bramnum][bramsize], UDTYPE operand2[bramnum][bramsize],
				   UDTYPE rp[N], UDTYPE srp[N], UDTYPE modulus_in, UDTYPE key_modulus_const_ratio1)
{
#pragma HLS INLINE off

	UDTYPE modulus = modulus_in;

	UDTYPE operanda[2][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	UDTYPE operandb[2][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	UDTYPE operanda_[2][corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	UDTYPE operandb_[2][corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	UDTYPE two_times_modulus = modulus << 1;
	UDTYPE one_u64 = 1;
	UDTYPE stepNum = bramsize;
	UDTYPE stageNum = STAGENUM;
	UDTYPE MEa;
	UDTYPE MEb;

	UDTYPE arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

	int stage1_max = STAGE1_MAX;

ntt_i_0to9:
	for (uint32_t i = 0; i < stage1_max; i++)
	{
		UDTYPE stepsize = bramsize >> (i + 1);

		uint32_t temp1 = 1 << i;
		uint32_t temp2 = stage1_max - i;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	ntt_i_0to9_j:
		for (uint32_t j = 0; j < stepNum; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand1 inter false
#pragma HLS DEPENDENCE variable = operand2 inter false

			UDTYPE w0 = rp[temp1 + (j >> temp2)];
			UDTYPE sw0 = srp[temp1 + (j >> temp2)];
			UDTYPE w00 = rp[temp1 + ((j + 1) >> temp2)];
			UDTYPE sw00 = srp[temp1 + ((j + 1) >> temp2)];

			for (int k = 0; k < corenum; k++)
			{
#pragma HLS UNROLL
				operanda[0][k] = operand1[k][MEa];
				operandb[0][k] = operand1[k][MEb];
				operanda[1][k] = operand2[k][MEa];
				operandb[1][k] = operand2[k][MEb];
				operanda_[0][k] = operand1[k + corenum][MEa];
				operandb_[0][k] = operand1[k + corenum][MEb];
				operanda_[1][k] = operand2[k + corenum][MEa];
				operandb_[1][k] = operand2[k + corenum][MEb];
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[0][m], operandb[0][m]);
				ntt_negacyclic_harvey_lazy_core(w0, sw0, modulus, two_times_modulus, operanda[1][m], operandb[1][m]);
			}
			j++;
			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				ntt_negacyclic_harvey_lazy_core(w00, sw00, modulus, two_times_modulus, operanda_[0][m], operandb_[0][m]);
				ntt_negacyclic_harvey_lazy_core(w00, sw00, modulus, two_times_modulus, operanda_[1][m], operandb_[1][m]);
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL

				operand1[m][MEa] = operanda[0][m];
				operand1[m][MEb] = operandb[0][m];
				operand2[m][MEa] = operanda[1][m];
				operand2[m][MEb] = operandb[1][m];
				operand1[m + corenum][MEa] = operanda_[0][m];
				operand1[m + corenum][MEb] = operandb_[0][m];
				operand2[m + corenum][MEa] = operanda_[1][m];
				operand2[m + corenum][MEb] = operandb_[1][m];
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}

	int Sindexa = 0;
	int Sindexb = 0;

ntt_i_8_11:
	for (uint32_t i = stage1_max; i < stageNum; i++)
	{

		int windex = 0;
		for (int k = 0; k < corenum; k++)
		{

			if (i == 8)
			{
				arrayWindex[k] = bramsize; //閿熸枻鎷烽敓鏂ゆ嫹閿熸枻鎷峰憭閿熸枻鎷烽敓鎻亷鎷烽敓锟�?
			}
			else if (i == 9)
			{
				arrayWindex[k] = (bramsize << 1) + (k >> 2);
			}
			else if (i == 10)
			{
				arrayWindex[k] = (bramsize << 2) + (k >> 1);
			}
			else
			{
				arrayWindex[k] = (bramsize << 3) + k;
			}
		}

		MEa = 0;
		MEb = 0;

		UDTYPE temp_par = 8 >> (11 - i);
	ntt_stage2_inner:
		for (uint32_t j = 0; j < stepNum; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand1 inter false

			UDTYPE temp_adder = temp_par * j;
			UDTYPE w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			UDTYPE sw[corenum];
#pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				int x = (arrayWindex[m] + temp_adder);

				//			w[m]=rp[x];
				//			sw[m]=srp[x];

				w[m] = rp[0];
				sw[m] = srp[0];
				windex++;
			}
			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operanda[0][m] = operand1[m][MEa];
				operandb[0][m] = operand1[m + corenum][MEb];

				operanda[1][m] = operand2[m][MEa];
				operandb[1][m] = operand2[m + corenum][MEb];
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				ntt_negacyclic_harvey_lazy_core(w[m], sw[m], modulus, two_times_modulus, operanda[0][m], operandb[0][m]);
				ntt_negacyclic_harvey_lazy_core(w[m], sw[m], modulus, two_times_modulus, operanda[1][m], operandb[1][m]);
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operand1[m][MEa] = operanda[0][m];
				operand1[m + corenum][MEb] = operandb[0][m];
				operand2[m][MEa] = operanda[1][m];
				operand2[m + corenum][MEb] = operandb[1][m];
			}

			MEa += 1;
			MEb += 1;
		}
	}

	UDTYPE half = modulus >> 1;
	for (size_t l = 0; l < bramsize; l++)
	{
#pragma HLS PIPELINE
		for (size_t r = 0; r < bramnum; r++)
		{
			operand1[r][l] = barrett_reduce_63(
				operand1[r][l] + half,
				modulus,
				key_modulus_const_ratio1);
			operand2[r][l] = barrett_reduce_63(
				operand2[r][l] + half,
				modulus,
				key_modulus_const_ratio1);
		}
	}
}
