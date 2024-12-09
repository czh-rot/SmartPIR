#include "top.h"

using namespace std;

UDTYPE add_uint64(
    UDTYPE operand1, UDTYPE operand2, UDTYPE* result)
{
    *result = operand1 + operand2;
    return UDTYPE(*result < operand1);   //如果*result 小于 operand1就输出1
}

UDTYPE add_uint64_generic(
    UDTYPE operand1, UDTYPE operand2, uint8_t carry,
    UDTYPE* result)
{
    operand1 += operand2;
    *result = operand1 + carry;
    return (operand1 < operand2) || (~operand1 < carry);
}

/*void add_poly_poly_coeffmod(const UDTYPE* operand1,
    const UDTYPE* operand2,
    const UDTYPE modulus, UDTYPE* result)
{

	UDTYPE coeff_count = N;
    const UDTYPE modulus_value = modulus;
    for (; coeff_count--; result++, operand1++, operand2++)
    {
#pragma HLS PIPELINE
        UDTYPE sum = *operand1 + *operand2;
        *result = sum - (modulus_value & UDTYPE(
            -DTYPE(sum >= modulus_value)));
    }
}*/

//void add_N(const UDTYPE operand1[BRAMNUM][BRAMSIZE],
//    const UDTYPE operand2[BRAMNUM][BRAMSIZE],
//    const UDTYPE modulus, UDTYPE result[BRAMNUM][BRAMSIZE])
//{
//    const UDTYPE modulus_value = modulus;
//    for (size_t i = 0; i < BRAMNUM; i++){
//    	for (size_t j = 0; j < BRAMSIZE;j++){
//#pragma HLS PIPELINE
//    		UDTYPE sum = operand1[i][j] + operand2[i][j];
//    		result[i][j] = sum - (modulus_value & UDTYPE(
//    		    -DTYPE(sum >= modulus_value)));
//    	}
//    }
//}


UDTYPE Hadd_single(UDTYPE operand1,UDTYPE operand2,UDTYPE modulus_value){
#pragma HLS INLINE off
	UDTYPE sum = operand1 + operand2;
	UDTYPE result;
	result = sum - (modulus_value & UDTYPE(
	    -DTYPE(sum >= modulus_value)));
	return result;
}

/*void Hadd_N (UDTYPE operand1[BRAMNUM][BRAMSIZE], UDTYPE operand2[BRAMNUM][BRAMSIZE], UDTYPE modulus)
{
#pragma HLS INLINE off
    const UDTYPE modulus_value = modulus;

	for (size_t j = 0; j < BRAMSIZE;j++){
		 for (size_t i = 0; i < BRAMNUM; i++){
#pragma HLS PIPELINE
			 operand1[i][j] = Hadd_single(operand1[i][j], operand2[i][j], modulus);
		 }
    }
}*/

uint8_t sub_uint64(
    UDTYPE operand1, UDTYPE operand2, UDTYPE* result)
{
#pragma HLS INLINE off
    *result = operand1 - operand2;
    return uint8_t(operand2 > operand1);
}

UDTYPE Hsub_single(UDTYPE operand1,UDTYPE operand2,UDTYPE modulus){
#pragma HLS INLINE off
	UDTYPE temp, result;
	DTYPE borrow = sub_uint64(operand1, operand2, &temp);
	result = temp + (modulus & UDTYPE(-borrow));
	return result;
}

uint8_t sub_uint64_generic(
    UDTYPE operand1, UDTYPE operand2,
    uint8_t borrow, UDTYPE* result)
{
    UDTYPE diff = operand1 - operand2;
    *result = diff - (borrow != 0);
    return (diff > operand1) || (diff < borrow);
}

UDTYPE sub_uint_uint_mod(
    UDTYPE operand1, UDTYPE operand2,
    UDTYPE modulus)
{
    UDTYPE temp;
    DTYPE borrow = sub_uint64_generic(operand1, operand2, 0, &temp);
    return UDTYPE(temp) +
        (modulus & UDTYPE(-borrow));
}

/*void sub_poly_poly_coeffmod(const UDTYPE* operand1,
    const UDTYPE* operand2, size_t coeff_count,
    const UDTYPE modulus, UDTYPE* result) {
    const UDTYPE modulus_value = modulus;
    for (; coeff_count--; result++, operand1++, operand2++)
    {
#pragma HLS PIPELINE
        UDTYPE temp_result;
        DTYPE borrow = sub_uint64(*operand1, *operand2, &temp_result);
        *result = temp_result + (modulus_value & UDTYPE(-borrow));
    }
}*/

/*void sub_8_1024(const UDTYPE operand1[BRAMNUM][BRAMSIZE],
    const UDTYPE operand2[BRAMNUM][BRAMSIZE],
    const UDTYPE modulus, UDTYPE result[BRAMNUM][BRAMSIZE])
{
	const UDTYPE modulus_value = modulus;
	for (size_t i = 0; i < 8; i++){
		for (size_t j = 0; j < BRAMSIZE;j++){
#pragma HLS PIPELINE
			UDTYPE temp_result;
			DTYPE borrow = sub_uint64(operand1[i][j], operand2[i][j], &temp_result);
			result[i][j] = temp_result + (modulus_value & UDTYPE(-borrow));
		}
	}
}*/



/*void multiply_uint64_hw64(
    UDTYPE operand1, UDTYPE operand2, UDTYPE* hw64)
{
    UDTYPE operand1_coeff_right = operand1 & 0x00000000FFFFFFFF;
    UDTYPE operand2_coeff_right = operand2 & 0x00000000FFFFFFFF;
    operand1 >>= 32;
    operand2 >>= 32;

    UDTYPE middle1 = operand1 * operand2_coeff_right;
    UDTYPE middle;
    UDTYPE left = operand1 * operand2 + (UDTYPE(add_uint64(
        middle1, operand2 * operand1_coeff_right, &middle)) << 32);
    UDTYPE right = operand1_coeff_right * operand2_coeff_right;
    UDTYPE temp_sum = (right >> 32) + (middle & 0x00000000FFFFFFFF);

    *hw64 = UDTYPE(
        left + (middle >> 32) + (temp_sum >> 32));
}*/

/*void multiply_uint64(
    UDTYPE operand1, UDTYPE operand2, UDTYPE result1BITWIDTH[2])
{

    UDTYPE operand1_coeff_right = operand1 & 0x00000000FFFFFFFF;
    UDTYPE operand2_coeff_right = operand2 & 0x00000000FFFFFFFF;
    operand1 >>= 32;
    operand2 >>= 32;

    UDTYPE middle1 = operand1 * operand2_coeff_right;
    UDTYPE middle;
    UDTYPE left = operand1 * operand2 + (UDTYPE(add_uint64(
        middle1, operand2 * operand1_coeff_right, &middle)) << 32);
    UDTYPE right = operand1_coeff_right * operand2_coeff_right;
    UDTYPE temp_sum = (right >> 32) + (middle & 0x00000000FFFFFFFF);

    result1BITWIDTH[1] = UDTYPE(
        left + (middle >> 32) + (temp_sum >> 32));
    result1BITWIDTH[0] = UDTYPE(
        (temp_sum << 32) | (right & 0x00000000FFFFFFFF));
}*/

/*UDTYPE barrett_reduce_63(
    UDTYPE input, const UDTYPE modulus, const UDTYPE modulus_ratio1)   //这个input其实可以63-bit的
{
//#pragma HLS ALLOCATION instances=mul limit=1 operation
#pragma HLS INLINE off
    UDTYPE tmp1[2];
    UDTYPE2 tmp2;
    tmp2 = input * modulus_ratio1;
    tmp1[0] = tmp2(BITWIDTH-1,0);
    tmp1[1] = tmp2(2*BITWIDTH-1,BITWIDTH);
    // Barrett subtraction
    tmp1[0] = input - tmp1[1] * modulus;

    return static_cast<UDTYPE>(tmp1[0]) -
        (modulus & static_cast<UDTYPE>(
            -static_cast<DTYPE>(tmp1[0] >= modulus)));

}*/

UDTYPE barrett_reduce_63(
    UDTYPE input, UDTYPE modulus, UDTYPE modulus_ratio1)   //这个input其实可以63-bit的
{
//#pragma HLS ALLOCATION instances=mul limit=1 operation
#pragma HLS INLINE off
    UDTYPE tmp1[2];
    UDTYPE2 tmp2;
    tmp2 = input * modulus_ratio1;
    tmp1[0] = tmp2(BITWIDTH-1,0);
    tmp1[1] = tmp2(2*BITWIDTH-1,BITWIDTH);
    // Barrett subtraction
    tmp1[0] = input - tmp1[1] * modulus;

    return (UDTYPE)(tmp1[0]) -
        (modulus & (UDTYPE)(
            -(UDTYPE)(tmp1[0] >= modulus)));

}


/*void barrett_reduce_N(UDTYPE operand1[BRAMNUM][BRAMSIZE],
	   UDTYPE modulus,UDTYPE modulus_ratio1){
#pragma HLS INLINE off

	for(int j=0;j<BRAMSIZE;j++){

			for(int i=0;i<BRAMNUM;i++){
#pragma HLS PIPELINE

			operand1[i][j]=barrett_reduce_63(
			    operand1[i][j], modulus,modulus_ratio1);

		}

	}

}*/
/*UDTYPE barrett_reduce_1BITWIDTH(
    UDTYPE* input, const UDTYPE modulus, UDTYPE modulus_ratio0, UDTYPE modulus_ratio1)
{

    // Reduces input using base 2^64 Barrett reduction
    // input allocation size must be 1BITWIDTH bits

    UDTYPE tmp1, tmp2[2], tmp3, carry;
    UDTYPE2 tmp4;

    // Multiply input and const_ratio
    // Round 1
    //multiply_uint64_hw64(input[0], modulus_ratio0, &carry);
    carry = (input[0] * modulus_ratio0)(2*BITWIDTH-1,BITWIDTH);

    tmp4 = input[0] * modulus_ratio1;
    tmp2[0] = tmp4(BITWIDTH-1,0);
    tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
    //multiply_uint64(input[0], modulus_ratio1, tmp2);
    tmp3 = tmp2[1] + add_uint64_generic(tmp2[0], carry, 0, &tmp1);

    // Round 2
    tmp4 = input[1] * modulus_ratio0;
    tmp2[0] = tmp4(BITWIDTH-1,0);
    tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
    //multiply_uint64(input[1], modulus_ratio0, tmp2);
    carry = tmp2[1] + add_uint64_generic(tmp1, tmp2[0], 0, &tmp1);

    // This is all we care about
    tmp1 = input[1] * modulus_ratio1 + tmp3 + carry;

    // Barrett subtraction
    tmp3 = input[0] - tmp1 * modulus;

    // One more subtraction is enough
    return UDTYPE(tmp3) -
        (modulus & UDTYPE(
            -DTYPE(tmp3 >= modulus)));
}*/

UDTYPE barrett_reduce_128(
    UDTYPE input0,UDTYPE input1, UDTYPE modulus, UDTYPE modulus_ratio0, UDTYPE modulus_ratio1)
{

    // Reduces input using base 2^64 Barrett reduction
    // input allocation size must be 1BITWIDTH bits

    UDTYPE tmp1, tmp2[2], tmp3, carry;
    UDTYPE2 tmp4;

    // Multiply input and const_ratio
    // Round 1
    //multiply_uint64_hw64(input[0], modulus_ratio0, &carry);
    carry = (input0 * modulus_ratio0)(2*BITWIDTH-1,BITWIDTH);

    tmp4 = input0 * modulus_ratio1;
    tmp2[0] = tmp4(BITWIDTH-1,0);
    tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
    //multiply_uint64(input[0], modulus_ratio1, tmp2);
    tmp3 = tmp2[1] + add_uint64_generic(tmp2[0], carry, 0, &tmp1);

    // Round 2
    tmp4 = input1 * modulus_ratio0;
    tmp2[0] = tmp4(BITWIDTH-1,0);
    tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
    //multiply_uint64(input[1], modulus_ratio0, tmp2);
    carry = tmp2[1] + add_uint64_generic(tmp1, tmp2[0], 0, &tmp1);

    // This is all we care about
    tmp1 = input1 * modulus_ratio1 + tmp3 + carry;

    // Barrett subtraction
    tmp3 = input0 - tmp1 * modulus;

    // One more subtraction is enough
    return UDTYPE(tmp3) -
        (modulus & UDTYPE(
            -DTYPE(tmp3 >= modulus)));
}

/*void multiply_poly_scalar_coeffmod(const UDTYPE* poly,
    size_t coeff_count, UDTYPE scalar, const UDTYPE modulus, const UDTYPE modulus_ratio0, const UDTYPE modulus_ratio1,
    UDTYPE* result)
{
    for (; coeff_count--; poly++, result++)
    {
#pragma HLS PIPELINE
        UDTYPE z[2], tmp1, tmp2[2], tmp3, carry;
        UDTYPE2 tmpz, tmp4;
		tmpz = (*poly) * (scalar);
		z[0] = tmpz(BITWIDTH-1,0);
		z[1] = tmpz(2*BITWIDTH-1,BITWIDTH);

        // Multiply input and const_ratio
        // Round 1
		carry = (z[0] * modulus_ratio0)(2*BITWIDTH-1,BITWIDTH);
		tmp4 = z[0] * modulus_ratio1;
		tmp2[0] = tmp4(BITWIDTH-1,0);
		tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
        tmp3 = tmp2[1] + add_uint64(tmp2[0], carry, &tmp1);

        // Round 2
        tmp4 = z[1] * modulus_ratio0;
		tmp2[0] = tmp4(BITWIDTH-1,0);
		tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
        carry = tmp2[1] + add_uint64(tmp1, tmp2[0], &tmp1);

        // This is all we care about
        tmp1 = z[1] * modulus_ratio1 + tmp3 + carry;

        // Barrett subtraction
        tmp3 = z[0] - tmp1 * modulus;

        // Claim: One more subtraction is enough
        *result = tmp3 - (modulus & UDTYPE(
            -DTYPE(tmp3 >= modulus)));
    }
}*/
/*void multiply_scalar_8_1024(const UDTYPE poly[BRAMNUM][BRAMSIZE],
    size_t coeff_count, UDTYPE scalar, const UDTYPE modulus, const UDTYPE modulus_ratio0, const UDTYPE modulus_ratio1,
    UDTYPE result[BRAMNUM][BRAMSIZE])
{
    for (size_t i = 0; i < 8; i++){
    	for (size_t j = 0; j < BRAMSIZE; j++){
#pragma HLS PIPELINE
        UDTYPE z[2], tmp1, tmp2[2], tmp3, carry;
        UDTYPE2 tmpz, tmp4;
		tmpz = poly[i][j] * scalar;
		z[0] = tmpz(BITWIDTH-1,0);
		z[1] = tmpz(2*BITWIDTH-1,BITWIDTH);

        // Multiply input and const_ratio
        // Round 1
		carry = (z[0] * modulus_ratio0)(2*BITWIDTH-1,BITWIDTH);
		tmp4 = z[0] * modulus_ratio1;
		tmp2[0] = tmp4(BITWIDTH-1,0);
		tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
        tmp3 = tmp2[1] + add_uint64(tmp2[0], carry, &tmp1);

        // Round 2
        tmp4 = z[1] * modulus_ratio0;
		tmp2[0] = tmp4(BITWIDTH-1,0);
		tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
        carry = tmp2[1] + add_uint64(tmp1, tmp2[0], &tmp1);

        // This is all we care about
        tmp1 = z[1] * modulus_ratio1 + tmp3 + carry;

        // Barrett subtraction
        tmp3 = z[0] - tmp1 * modulus;

        // Claim: One more subtraction is enough
        result[i][j] = tmp3 - (modulus & UDTYPE(
            -DTYPE(tmp3 >= modulus)));

    	}
    }
}*/

/*UDTYPE multiply_coeffmod(UDTYPE operand1,UDTYPE operand2,
    const UDTYPE modulus, const UDTYPE modulus_ratio0, const UDTYPE modulus_ratio1)
{
#pragma HLS INLINE off
	UDTYPE result = 0;
	UDTYPE z[2], tmp1, tmp2[2], tmp3, carry;
	UDTYPE2 tmpz, tmp4;
	tmpz = (operand1) * (operand2);
	z[0] = tmpz(BITWIDTH-1,0);
	z[1] = tmpz(2*BITWIDTH-1,BITWIDTH);

	// Multiply input and const_ratio
	// Round 1
	carry = (z[0] * modulus_ratio0)(2*BITWIDTH-1,BITWIDTH);
	tmp4 = z[0] * modulus_ratio1;
	tmp2[0] = tmp4(BITWIDTH-1,0);
	tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
	tmp3 = tmp2[1] + add_uint64(tmp2[0], carry, &tmp1);

	// Round 2
	tmp4 = z[1] * modulus_ratio0;
	tmp2[0] = tmp4(BITWIDTH-1,0);
	tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
	carry = tmp2[1] + add_uint64(tmp1, tmp2[0], &tmp1);

	// This is all we care about
	tmp1 = z[1] * modulus_ratio1 + tmp3 + carry;

	// Barrett subtraction
	tmp3 = z[0] - tmp1 * modulus;

	// Claim: One more subtraction is enough
	result = tmp3 - (modulus & UDTYPE(
		-DTYPE(tmp3 >= modulus)));

	return result;
}*/



/*UDTYPE HMult_single(UDTYPE operand1,UDTYPE operand2,
    UDTYPE modulus,UDTYPE modulus_ratio[2])
{
#pragma HLS INLINE off

	UDTYPE z[2], tmp1, tmp2[2], tmp3, carry;
	UDTYPE2 tmpz, tmp4;
	tmpz = (operand1) * (operand2);
	z[0] = tmpz(BITWIDTH-1,0);
	z[1] = tmpz(2*BITWIDTH-1,BITWIDTH);

	// Multiply input and const_ratio
	// Round 1
	carry = (z[0] * modulus_ratio[0])(2*BITWIDTH-1,BITWIDTH);
	tmp4 = z[0] * modulus_ratio[1];
	tmp2[0] = tmp4(BITWIDTH-1,0);
	tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
	tmp3 = tmp2[1] + add_uint64(tmp2[0], carry, &tmp1);

	// Round 2
	tmp4 = z[1] * modulus_ratio[0];
	tmp2[0] = tmp4(BITWIDTH-1,0);
	tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
	carry = tmp2[1] + add_uint64(tmp1, tmp2[0], &tmp1);

	// This is all we care about
	tmp1 = z[1] * modulus_ratio[1] + tmp3 + carry;

	// Barrett subtraction
	tmp3 = z[0] - tmp1 * modulus;

	// Claim: One more subtraction is enough
	UDTYPE result = tmp3 - (modulus & UDTYPE(
		-DTYPE(tmp3 >= modulus)));

	return result;

}*/

UDTYPE HMult_single_return(UDTYPE operand1,UDTYPE operand2,
    UDTYPE modulus,UDTYPE modulus_ratio_0, UDTYPE modulus_ratio_1)
{
//#pragma HLS ALLOCATION instances=mul limit=6 operation
#pragma HLS INLINE off

	UDTYPE z[2], tmp1, tmp2[2], tmp3, carry;
	UDTYPE2 tmpz, tmp4;
	tmpz = (operand1) * (operand2);
	z[0] = tmpz(BITWIDTH-1,0);
	z[1] = tmpz(2*BITWIDTH-1,BITWIDTH);

	// Multiply input and const_ratio
	// Round 1
	carry = (z[0] * modulus_ratio_0)(2*BITWIDTH-1,BITWIDTH);
	tmp4 = z[0] * modulus_ratio_1;
	tmp2[0] = tmp4(BITWIDTH-1,0);
	tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
	tmp3 = tmp2[1] + add_uint64(tmp2[0], carry, &tmp1);

	// Round 2
	tmp4 = z[1] * modulus_ratio_0;
	tmp2[0] = tmp4(BITWIDTH-1,0);
	tmp2[1] = tmp4(2*BITWIDTH-1,BITWIDTH);
	carry = tmp2[1] + add_uint64(tmp1, tmp2[0], &tmp1);

	// This is all we care about
	tmp1 = z[1] * modulus_ratio_1 + tmp3 + carry;

	// Barrett subtraction
	tmp3 = z[0] - tmp1 * modulus;

	// Claim: One more subtraction is enough
	UDTYPE result = tmp3 - (modulus & UDTYPE(
		-DTYPE(tmp3 >= modulus)));

	return result;

}

