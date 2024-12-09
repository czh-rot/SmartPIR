#include "top.h"
using namespace std;

//NTT_core
void ntt_negacyclic_harvey_lazy_core(UDTYPE W, UDTYPE Wprime, UDTYPE modulus, UDTYPE two_times_modulus,
		UDTYPE &operanda, UDTYPE &operandb){
#pragma HLS INLINE off
	UDTYPE currX = operanda - (two_times_modulus & (UDTYPE)(-(DTYPE)(operanda >= two_times_modulus)));
	UDTYPE2 result = Wprime * operandb;
	UDTYPE Q = result(2*BITWIDTH-1, BITWIDTH);
	Q = W * (operandb) - Q * modulus;
	Q = Q - (modulus & UDTYPE(-DTYPE(Q >= modulus)));
	operanda = currX + Q;
	operandb = currX + (two_times_modulus - Q);
}

void ntt_core(UDTYPE W, UDTYPE Wprime, UDTYPE modulus, UDTYPE two_times_modulus,
		UDTYPE operanda_in, UDTYPE operandb_in, UDTYPE &operanda, UDTYPE &operandb){
#pragma HLS INLINE off
//#pragma HLS ALLOCATION instances=Mul limit=2 core
	UDTYPE currX, Q;
//#pragma HLS RESOURCE variable=currX core=AddSub_DSP
//#pragma HLS RESOURCE variable=Q core=AddSub_DSP
	UDTYPE operanda_temp, operandb_temp;

	operanda_temp = operanda_in;
	operandb_temp = operandb_in;

	currX = operanda_temp - (two_times_modulus & (UDTYPE)(-(DTYPE)(operanda_temp >= two_times_modulus)));
	UDTYPE2 result = Wprime * operandb_temp;
	Q = result(2*BITWIDTH-1, BITWIDTH);
	Q = W * (operandb_temp) - Q * modulus;
	Q = Q - (modulus & UDTYPE(-DTYPE(Q >= modulus)));
	operanda = currX + Q;
	operandb = currX + (two_times_modulus - Q);
}



