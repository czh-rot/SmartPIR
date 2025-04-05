// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Alişah Özcan

#include "heongpu.cuh"

#include <string>
#include <iomanip>
#include <omp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

int main(int argc, char* argv[])
{
    cudaSetDevice(0);

    std::size_t poly_modulus_degrees = 8192;
    std::size_t plain_modulus = 65537;

    heongpu::Parameters context(
        heongpu::scheme_type::bfv,
        heongpu::keyswitching_type::KEYSWITCHING_METHOD_I);
    context.set_poly_modulus_degree(poly_modulus_degrees);
    context.set_default_coeff_modulus(1);
    context.set_plain_modulus(plain_modulus);
    context.generate();

    heongpu::HEKeyGenerator keygen(context);
    heongpu::Secretkey secret_key(context);
    keygen.generate_secret_key(secret_key);

    heongpu::Publickey public_key(context);
    keygen.generate_public_key(public_key, secret_key);

    heongpu::Relinkey relin_key(context);
    keygen.generate_relin_key(relin_key, secret_key);

    std::vector<int> custom_key_index = {1};
    heongpu::Galoiskey galois_key(context, custom_key_index);
    keygen.generate_galois_key(galois_key, secret_key);

    heongpu::HEEncoder encoder(context);
    heongpu::HEEncryptor encryptor(context, public_key);
    heongpu::HEDecryptor decryptor(context, secret_key);
    heongpu::HEArithmeticOperator operators(context, encoder);

    // Easily, e assume the database contain N key-value items
    // assume each entry is 64 bits, num = 8192.
    
    heongpu::HostVector<heongpu::Plaintext> Value(4);
    
    for (int i = 0; i < 4; i++) {
        heongpu::HostVector<uint64_t> messagev(poly_modulus_degrees, 1);
        encoder.encode(Value[i], messagev);
    }
    
    heongpu::HostVector<uint64_t> query(poly_modulus_degrees, 1);
    query[0] = 1;
    heongpu::HostVector<heongpu::Plaintext> pQ;
    heongpu::HostVector<heongpu::Plaintext> cQ;
    encoder.encode(pQ,query);
        encryptor.encrypt(cQ, pQ);

    std::cout << "INFO: The server use the I to extract the target value." << std::endl;
    heongpu::HostVector<heongpu::Ciphertext> cV(4);
    for (int i = 0; i < 4; i++) {
        operators.multiply_plain(cQ, Value[i], cV[i]);
    }
    heongpu::Ciphertext Ans(context);
    for (int i = 0; i < t; i++) {
        operators.rotate_rows(cV[i], cV[i], galois_key, 1);
        operators.add(Ans, cV[i], Ans);
    }
    
    return EXIT_SUCCESS;
}