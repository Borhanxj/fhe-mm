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

    // ---------------- Generate two matrices ----------------
    cout << "Generating random matrices A and B" << endl;
    size_t n = 16;
    auto A = utils::random_matrix(n, n);
    auto B = utils::random_matrix(n, n);

    // ---------------- Pack and encrypt A diagonally ----------------
    cout << "Encrypting diagonals of A (with replication)" << endl;
    vector<Ciphertext> enc_matrix_A;
    enc_matrix_A.reserve(n);
    for (size_t d = 0; d < n; ++d)
    {
        vector<double> diagonal(slots, 0.0);
        for (size_t i = 0; i < n; ++i) {
            diagonal[i] = A[i][(i + d) % n];
        }
        vector<double> replicated_diag(slots);
        for(size_t i = 0; i < slots; ++i) {
            replicated_diag[i] = diagonal[i % n];
        }
        Plaintext plain_diag;
        encoder.encode(replicated_diag, scale, plain_diag);
        Ciphertext enc_diag;
        encryptor.encrypt(plain_diag, enc_diag);
        enc_matrix_A.push_back(enc_diag);
    }

    // ---------------- Encrypt B column-wise ----------------
    cout << "Encrypting columns of B (with replication)" << endl;
    vector<Ciphertext> enc_matrix_B;
    enc_matrix_B.reserve(n);
    for (size_t c = 0; c < n; ++c)
    {
        vector<double> col(n);
        for (size_t r = 0; r < n; ++r)
            col[r] = B[r][c];

        vector<double> replicated_col(slots);
        for (size_t i = 0; i < slots; ++i) {
            replicated_col[i] = col[i % n];
        }
        Plaintext plain_col;
        encoder.encode(replicated_col, scale, plain_col);
        Ciphertext enc_col;
        encryptor.encrypt(plain_col, enc_col);
        enc_matrix_B.push_back(enc_col);
    }

    // ---------------- Run encrypted mat × mat ----------------
    cout << "Running diagonal_encmat_encmat (encrypted × encrypted)" << endl;
    vector<Ciphertext> encrypted_result =
        matvec::diagonal_encmat_encmat(enc_matrix_A, enc_matrix_B, n,
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
        {
            result_matrix[r][c] = decoded[r];
        }
    }
    
    // ---------------- Automated Checking ----------------
    cout << endl << "Verifying correctness of the result..." << endl;
    auto expected = utils::expected_result(A, B);
    double tolerance = 0.01;
    utils::check_matrix_approximation(result_matrix, expected, tolerance);

    // ---------------- Profiling ----------------
    cout << endl << "Writing profiling data..." << endl;
    print_profiling("../benchmarks/results/diagonal_encmat_encmat.txt");

    cout << "Done!" << endl;
    return 0;
}
