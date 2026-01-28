#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "seal/seal.h"
#include "utils.hpp"
#include "matvec.hpp"
#include "profiler.hpp"

using namespace std;
using namespace seal;

int main() {
    // ---------------- Setup ----------------
    cout << "Setting up the environment" << endl;
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(8192);
    parms.set_coeff_modulus(CoeffModulus::Create(8192, {50, 40, 40, 40, 40}));
    double scale = pow(2.0, 40);
    SEALContext context(parms);

    // ---------------- Keys ----------------
    cout << "Generating keys" << endl;
    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    GaloisKeys galois_keys;
    keygen.create_galois_keys(galois_keys);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    // ---------------- Helpers ----------------
    cout << "Initializing helpers" << endl;
    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);
    CKKSEncoder encoder(context);
    size_t slots = encoder.slot_count();

    // ---------------- Generate matrix and vector ----------------
    cout << "Generating random matrix and vector" << endl;
    size_t rows = 16;
    size_t cols = 16;
    auto matrix = utils::random_matrix(rows, cols);
    auto vec = utils::random_vector(cols);
    vector<double> expected = utils::expected_result(matrix, vec);
    double tolerance = 0.01;

    // ---------------- Encode and Encrypt Vector ----------------
    cout << "Encoding and Encrypting vector (with replication)" << endl;
    vector<double> replicated_vec(slots);
    for(size_t i = 0; i < slots; ++i) {
        replicated_vec[i] = vec[i % cols];
    }
    Plaintext plain_vec;
    encoder.encode(replicated_vec, scale, plain_vec);
    Ciphertext enc_vec;
    encryptor.encrypt(plain_vec, enc_vec);

    //================================================================
    // TEST 1: Plaintext Matrix × Encrypted Vector
    //================================================================
    cout << "\n--- Running Test 1: Plaintext Matrix x Encrypted Vector ---" << endl;
    Ciphertext plain_enc_result = 
        matvec::diagonal_matvec(matrix, enc_vec, cols, evaluator, encoder, galois_keys, scale);

    cout << "Decrypting and decoding result..." << endl;
    Plaintext plain_result_pt;
    decryptor.decrypt(plain_enc_result, plain_result_pt);
    vector<double> decoded_plain_enc;
    encoder.decode(plain_result_pt, decoded_plain_enc);

    cout << "Verifying correctness..." << endl;
    utils::check_approximation(decoded_plain_enc, expected, tolerance);
    
    cout << "Writing profiling data for Test 1..." << endl;
    print_profiling("../benchmarks/results/diagonal_plain_matvec.txt");
    clear_profiling();

    //================================================================
    // TEST 2: Encrypted Matrix × Encrypted Vector
    //================================================================
    cout << "\n--- Running Test 2: Encrypted Matrix x Encrypted Vector ---" << endl;
    
    cout << "Encrypting matrix diagonals (with replication)..." << endl;
    vector<Ciphertext> enc_matrix_diags;
    enc_matrix_diags.reserve(cols);
    for (size_t d = 0; d < cols; ++d) {
        vector<double> diagonal(slots, 0.0);
        // Extract the d-th diagonal of the matrix
        for (size_t i = 0; i < cols; ++i) {
            diagonal[i] = matrix[i][(i + d) % cols];
        }
        // Replicate the diagonal to fill the slots for encoding
        vector<double> replicated_diag(slots);
        for(size_t i = 0; i < slots; ++i) {
            replicated_diag[i] = diagonal[i % cols];
        }
        Plaintext plain_diag;
        encoder.encode(replicated_diag, scale, plain_diag);
        Ciphertext enc_diag;
        encryptor.encrypt(plain_diag, enc_diag);
        enc_matrix_diags.push_back(enc_diag);
    }

    Ciphertext enc_enc_result = 
        matvec::diagonal_encmat_encvec(enc_matrix_diags, enc_vec, cols, evaluator, encoder, galois_keys, relin_keys, scale);

    cout << "Decrypting and decoding result..." << endl;
    Plaintext enc_result_pt;
    decryptor.decrypt(enc_enc_result, enc_result_pt);
    vector<double> decoded_enc_enc;
    encoder.decode(enc_result_pt, decoded_enc_enc);

    cout << "Verifying correctness..." << endl;
    utils::check_approximation(decoded_enc_enc, expected, tolerance);
    
    cout << "Writing profiling data for Test 2..." << endl;
    print_profiling("../benchmarks/results/diagonal_enc_matvec.txt");
    
    cout << "\nDone!" << endl;
    return 0;
}
