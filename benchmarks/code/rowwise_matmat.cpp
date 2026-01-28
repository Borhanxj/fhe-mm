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

int main()
{
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

    // ---------------- Generate matrices ----------------
    cout << "Generating random matrices A and B" << endl;
    size_t n = 16; 
    auto A = utils::random_matrix(n, n);
    auto B = utils::random_matrix(n, n);

    // ---------------- Encrypt A row-wise ----------------
    cout << "Encrypting rows of A (row-wise replication)" << endl;
    vector<Ciphertext> enc_matrix_A;
    enc_matrix_A.reserve(n);
    for (size_t r = 0; r < n; ++r)
    {
        vector<double> replicated_row(slots);
        for (size_t j = 0; j < slots; ++j) {
            replicated_row[j] = A[r][j % n];
        }

        Plaintext plain_row;
        encoder.encode(replicated_row, scale, plain_row);
        Ciphertext enc_row;
        encryptor.encrypt(plain_row, enc_row);
        enc_matrix_A.push_back(enc_row);
    }

    // ---------------- Encrypt B column-wise ----------------
    cout << "Encrypting columns of B (column-wise replication)" << endl;
    vector<Ciphertext> enc_matrix_B;
    enc_matrix_B.reserve(n);
    for (size_t c = 0; c < n; ++c)
    {
        vector<double> replicated_col(slots);
        for (size_t i = 0; i < slots; ++i) {
            replicated_col[i] = B[i % n][c];
        }

        Plaintext plain_col;
        encoder.encode(replicated_col, scale, plain_col);
        Ciphertext enc_col;
        encryptor.encrypt(plain_col, enc_col);
        enc_matrix_B.push_back(enc_col);
    }

    // ---------------- Run encrypted mat × mat ----------------
    cout << "Running rowwise_encmat_encmat (encrypted × encrypted)" << endl;
    vector<Ciphertext> encrypted_result =
        matvec::rowwise_encmat_encmat(enc_matrix_A, enc_matrix_B, n,
                                      evaluator, encoder, galois_keys, relin_keys, scale);

    // ---------------- Decrypt and decode ----------------
    cout << "Decrypting result..." << endl;
    vector<vector<double>> result_matrix(n, vector<double>(n));
    for (size_t c = 0; c < n; ++c)
    {
        Plaintext plain_res;
        decryptor.decrypt(encrypted_result[c], plain_res);
        vector<double> decoded;
        encoder.decode(plain_res, decoded);

        for (size_t r = 0; r < n; ++r)
            result_matrix[r][c] = decoded[r];
    }

    // ---------------- Automated Checking ----------------
    cout << endl << "Verifying correctness of the result..." << endl;
    auto expected = utils::expected_result(A, B);
    double tolerance = 0.01;
    utils::check_matrix_approximation(result_matrix, expected, tolerance);

    // ---------------- Profiling ----------------
    cout << endl << "Writing profiling data..." << endl;
    print_profiling("../benchmarks/results/rowwise_encmat_encmat.txt");

    cout << "Done!" << endl;
    return 0;
}
