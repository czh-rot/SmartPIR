// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "examples.h"
#include "seal/seal.h"
#include "seal/util/polyarithsmallmod.h"
#include <chrono>
#include <bitset>
#include <boost/multiprecision/cpp_int.hpp>
#include <gmpxx.h>


using namespace std;
using namespace seal;
using namespace seal::util;
namespace mp = boost::multiprecision;

// The workflow of the SmartPIR example is as follows:
void example_bfv_basics()
{   

    print_example_banner("Example: BFV Basics");
    EncryptionParameters parms(scheme_type::bfv);
    size_t N = 32768; // BFV degree
    size_t poly_modulus_degree = 32768;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 17)); // 65537
    SEALContext context(parms);
    print_line(__LINE__);
    print_parameters(context);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    BatchEncoder encoder_(context);
    Decryptor decryptor(context, secret_key);

    std::random_device rd;  
    std::mt19937_64 gen(rd()); 
    std::uniform_int_distribution<uint64_t> dist;

    cout << "INFO: Input the number of items (e.g., 2^15 (32768)):" << endl;
    size_t num;
    cin >> num;
    cout << "INFO: Input the bit width (e.g., 64):" << endl;
    size_t l;
    cin >> l;
    cout << "INFO: Input the query id (e.g., 1, in [0,num-1]):" << endl;
    size_t id;
    cin >> id;

    // To simply, we generate N key-value pairs, and each key and value is 64 bits with random value
    cout << "INFO: we generate " << num <<" key-value pairs, and each key and value is 64 bits." << endl;
    size_t dim1 = num/N;
    cout << "INFO: They are divided into " << dim1 << "(i.e., num/N) partitions." << endl;
    vector<vector<vector<uint64_t >>> value(dim1, vector<vector<uint64_t >>(l/16, vector<uint64_t >(N)));


    for (uint64_t  i = 0; i < dim1; ++i) {
      for (uint64_t  j = 0; j < l/16; ++j) {
        for (uint64_t  k = 0; k < N; ++k) {
            cvalue[i][j][k] = dist(gen);
        }
      }
    }

    cout << "INFO: The server then encodes all values into plaintexts." << endl;

    vector<vector<Plaintext>> PValue(dim1, vector<Plaintext>(l/16));

    // Encode Key-Value into plaintexts
    for (uint64_t  i = 0; i < dim1; i++) {
      Plaintext c;
        for (uint64_t  j = 0; j < l/16; j++) {
            encoder_.encode(value[i][j], PValue[i][j]);
        }
    }

    // query, asumme query 1 st
    vector<Ciphertext> cquery(num/N);
    vector<Plaintext> pquery(num/N);
    vector<vector<uint64_t>> queries(num/N, vector<uint64_t>(N, 0));
    int idxVec = id / N;
    int idxEle = id % N;
    queries[idxVec][idxEle] = 1;
    for (int i = 0; i < num/N; i++) {
        encoder_.encode(queries[i], pquery[i]);
        encryptor.encrypt(pquery[i], cquery[i]);
    }

    cout << "INFO: The server use the I to extract the target value." << endl;
    vector<vector<Ciphertext>> CValue(dim1, vector<Ciphertext>(4));
    for (uint64_t  i = 0; i < dim1; i++) {
      for (uint64_t  j = 0; j < l/16; j++) {
        evaluator.multiply_plain(cquery[i], PValue[i][j], CValue[i][j]);
      }
    }

    // cout << "INFO: Noise budget of CValue: " << decryptor.invariant_noise_budget(CValue[0][0]) << endl;

    vector<vector<Ciphertext>> CValueTransposed(4, vector<Ciphertext>(dim1));
    // T CValue
    for (uint64_t  i = 0; i < dim1; i++) {
        for (uint64_t  j = 0; j < 4; j++) {
            CValueTransposed[j][i] = CValue[i][j];
        }
    }

    vector<Ciphertext> Value(4);
    for (uint64_t  j = 0; j < 4; j++) {
        evaluator.add_many(CValueTransposed[j], Value[j]);  // Ensure Value[j] is a vector
    }

    for (uint64_t  i = 1; i < 4; i++) {
      evaluator.rotate_rows_inplace(Value[i], i, gal_keys);
    }

    // cout << "INFO: Noise budget of Value: " << decryptor.invariant_noise_budget(Value[0]) << endl;

    Ciphertext Ans;
    evaluator.add_many(Value, Ans);

    Plaintext ans;
    decryptor.decrypt(Ans, ans);
    vector<uint64_t > ans_plain;
    encoder_.decode(ans, ans_plain);
    
    cout << "INFO: Noise budget of Ans: " << decryptor.invariant_noise_budget(Ans) << endl;
    cout << "INFO: The server obtain the Ans and send it back to the client." << endl;
    cout << "INFO: A private query is finished." << endl;
    return;


}
