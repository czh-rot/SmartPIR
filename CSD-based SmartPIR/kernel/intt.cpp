#include "top.h"
using namespace std;

void inverse_ntt_negacyclic_harvey_lazy_core(UDTYPE W, UDTYPE Wprime, UDTYPE modulus, UDTYPE two_times_modulus,
		UDTYPE &operanda, UDTYPE &operandb){
#pragma HLS INLINE off
	UDTYPE T= two_times_modulus - operandb + operanda;
	UDTYPE currX = operanda+operandb - (two_times_modulus & (UDTYPE)(-(DTYPE)((operanda<<1) >= T)));
	operanda=(currX + (modulus &(UDTYPE)(-(UDTYPE)(T & 1)))) >> 1;
	UDTYPE2 result = Wprime * T;
	UDTYPE Q = result(2*BITWIDTH-1, BITWIDTH);
	operandb=W * T - Q * modulus;

}

// inverse_ntt_negacyclic_harvey_lazy_core
void inverse_ntt_core(UDTYPE W, UDTYPE Wprime, UDTYPE modulus, UDTYPE two_times_modulus,
		UDTYPE operanda_in, UDTYPE operandb_in, UDTYPE &operanda, UDTYPE &operandb){
//#pragma HLS ALLOCATION instances=Mul limit=2 core
#pragma HLS INLINE off
	UDTYPE T, currX;
	UDTYPE operanda_temp, operandb_temp;
//#pragma HLS RESOURCE variable=T core=AddSub_DSP
//#pragma HLS RESOURCE variable=currX core=AddSub_DSP

	operanda_temp = operanda_in;
	operandb_temp = operandb_in;

	T= two_times_modulus - operandb_temp + operanda_temp;
//	cout << "modulus = " << modulus << endl;
//	cout << "*U = " << operanda_temp << endl;
//	cout << "*V = " << operandb_temp << endl;
//	cout << "two_times_modulus = " << two_times_modulus << endl;
//	cout << "T = " << T << endl;

	currX = operanda_temp + operandb_temp - (two_times_modulus & (UDTYPE)(-(DTYPE)((operanda_temp<<1) >= T)));
	//cout << "currU = " << currX << endl;

	operanda = (currX + (modulus &(UDTYPE)(-(UDTYPE)(T & 1)))) >> 1;
	UDTYPE2 result = Wprime * T;
	UDTYPE Q = result(2*BITWIDTH-1, BITWIDTH);

	//cout << "H = " << Q << endl;
	operandb=W * T - Q * modulus;

	operandb = operandb - (modulus & UDTYPE(-DTYPE(operandb >= modulus)));

}
