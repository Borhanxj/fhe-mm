#ifndef MATVEC_HPP
#define MATVEC_HPP

#include "profiler.hpp"
#include "seal/seal.h"
#include <vector>

namespace matvec {

using namespace std;
using namespace seal;

//================================================================================
//                               DECLARATIONS
//================================================================================

Ciphertext rowwise_matvec(const vector<vector<double>> &plain_matrix,
                          const Ciphertext &enc_vector, const size_t dimension,
                          Evaluator &evaluator, CKKSEncoder &encoder,
                          GaloisKeys &galois_keys, double scale);

Ciphertext diagonal_matvec(const vector<vector<double>> &plain_matrix,
                           const Ciphertext &enc_vector, const size_t dimension,
                           Evaluator &evaluator, CKKSEncoder &encoder,
                           GaloisKeys &galois_keys, double scale);

Ciphertext rowwise_encmat_encvec(const vector<Ciphertext> &enc_matrix_A,
                                 const Ciphertext &enc_vector_b,
                                 const size_t dimension, Evaluator &evaluator,
                                 CKKSEncoder &encoder, GaloisKeys &galois_keys,
                                 RelinKeys &relin_keys, double scale);

vector<Ciphertext> rowwise_encmat_encmat(const vector<Ciphertext> &enc_matrix_A,
                                         const vector<Ciphertext> &enc_matrix_B,
                                         const size_t dimension,
                                         Evaluator &evaluator,
                                         CKKSEncoder &encoder,
                                         GaloisKeys &galois_keys,
                                         RelinKeys &relin_keys, double scale);

Ciphertext diagonal_encmat_encvec(const vector<Ciphertext> &enc_matrix_A,
                                  const Ciphertext &enc_vector_x,
                                  const size_t dimension, Evaluator &evaluator,
                                  CKKSEncoder &encoder, GaloisKeys &galois_keys,
                                  RelinKeys &relin_keys, double scale);

vector<Ciphertext>
diagonal_encmat_encmat(const vector<Ciphertext> &enc_matrix_A,
                       const vector<Ciphertext> &enc_matrix_B,
                       const size_t dimension, Evaluator &evaluator,
                       CKKSEncoder &encoder, GaloisKeys &galois_keys,
                       RelinKeys &relin_keys, double scale);

//================================================================================
//                              IMPLEMENTATIONS
//================================================================================

Ciphertext rowwise_matvec(const vector<vector<double>> &plain_matrix,
                          const Ciphertext &enc_vector, const size_t dimension,
                          Evaluator &evaluator, CKKSEncoder &encoder,
                          GaloisKeys &galois_keys, double scale) {
  Ciphertext result;
  size_t n_rows = plain_matrix.size();
  size_t slots = encoder.slot_count();
  bool first = true;

  for (size_t row = 0; row < n_rows; ++row) {
    Plaintext plain_row;
    PROFILE_START("Encode Plain Row");
    encoder.encode(plain_matrix[row], scale, plain_row);
    PROFILE_END("Encode Plain Row");

    Ciphertext temp;
    PROFILE_START("Multiply Plain");
    evaluator.multiply_plain(enc_vector, plain_row, temp);
    evaluator.rescale_to_next_inplace(temp);
    PROFILE_END("Multiply Plain");

    for (size_t step = 1; step < dimension; step <<= 1) {
      Ciphertext rotated;
      PROFILE_START("Rotate+Add (Sum)");
      evaluator.rotate_vector(temp, step, galois_keys, rotated);
      evaluator.add_inplace(temp, rotated);
      PROFILE_END("Rotate+Add (Sum)");
    }

    vector<double> mask(slots, 0.0);
    mask[row] = 1.0;
    Plaintext mask_plain;
    PROFILE_START("Encode Mask");
    encoder.encode(mask, scale, mask_plain);
    PROFILE_END("Encode Mask");

    PROFILE_START("Masking Mult");
    evaluator.mod_switch_to_inplace(mask_plain, temp.parms_id());
    evaluator.multiply_plain_inplace(temp, mask_plain);
    evaluator.rescale_to_next_inplace(temp);
    PROFILE_END("Masking Mult");

    if (first) {
      result = temp;
      first = false;
    } else {
      PROFILE_START("Result Add");
      evaluator.add_inplace(result, temp);
      PROFILE_END("Result Add");
    }
  }
  return result;
}

Ciphertext diagonal_matvec(const vector<vector<double>> &plain_matrix,
                           const Ciphertext &enc_vector, const size_t dimension,
                           Evaluator &evaluator, CKKSEncoder &encoder,
                           GaloisKeys &galois_keys, double scale) {
  Ciphertext result;
  size_t slots = encoder.slot_count();
  bool first = true;

  for (size_t d = 0; d < dimension; ++d) {
    vector<double> diagonal(slots, 0.0);
    for (size_t i = 0; i < dimension; ++i) {
      size_t j = (i + d) % dimension;
      diagonal[i] = plain_matrix[i][j];
    }

    Plaintext plain_diag;
    PROFILE_START("Encode Diagonal");
    encoder.encode(diagonal, scale, plain_diag);
    PROFILE_END("Encode Diagonal");

    Ciphertext rotated_vec;
    if (d > 0) {
      PROFILE_START("Rotate Vector");
      evaluator.rotate_vector(enc_vector, static_cast<int>(d), galois_keys,
                              rotated_vec);
      PROFILE_END("Rotate Vector");
    } else {
      rotated_vec = enc_vector;
    }

    Ciphertext prod;
    PROFILE_START("Multiply Plain");
    evaluator.multiply_plain(rotated_vec, plain_diag, prod);
    evaluator.rescale_to_next_inplace(prod);
    PROFILE_END("Multiply Plain");

    if (first) {
      result = prod;
      first = false;
    } else {
      PROFILE_START("Result Add");
      evaluator.add_inplace(result, prod);
      PROFILE_END("Result Add");
    }
  }
  return result;
}

Ciphertext rowwise_encmat_encvec(const vector<Ciphertext> &enc_matrix_A,
                                 const Ciphertext &enc_vector_b,
                                 const size_t dimension, Evaluator &evaluator,
                                 CKKSEncoder &encoder, GaloisKeys &galois_keys,
                                 RelinKeys &relin_keys, double scale) {
  size_t n_rows = enc_matrix_A.size();
  size_t slots = encoder.slot_count();
  Ciphertext result;
  bool first = true;

  for (size_t i = 0; i < n_rows; ++i) {
    Ciphertext row_result;
    PROFILE_START("Homomorphic Mult");
    evaluator.multiply(enc_matrix_A[i], enc_vector_b, row_result);
    evaluator.relinearize_inplace(row_result, relin_keys);
    evaluator.rescale_to_next_inplace(row_result);
    PROFILE_END("Homomorphic Mult");

    for (size_t step = 1; step < dimension; step <<= 1) {
      Ciphertext rotated;
      PROFILE_START("Rotate+Add (Sum)");
      evaluator.rotate_vector(row_result, step, galois_keys, rotated);
      evaluator.add_inplace(row_result, rotated);
      PROFILE_END("Rotate+Add (Sum)");
    }

    vector<double> mask(slots, 0.0);
    mask[i] = 1.0;
    Plaintext mask_plain;
    PROFILE_START("Encode Mask");
    encoder.encode(mask, scale, mask_plain);
    PROFILE_END("Encode Mask");

    PROFILE_START("Masking Mult");
    evaluator.mod_switch_to_inplace(mask_plain, row_result.parms_id());
    evaluator.multiply_plain_inplace(row_result, mask_plain);
    evaluator.rescale_to_next_inplace(row_result);
    PROFILE_END("Masking Mult");

    if (first) {
      result = row_result;
      first = false;
    } else {
      PROFILE_START("Result Add");
      evaluator.add_inplace(result, row_result);
      PROFILE_END("Result Add");
    }
  }
  return result;
}

vector<Ciphertext> rowwise_encmat_encmat(const vector<Ciphertext> &enc_matrix_A,
                                         const vector<Ciphertext> &enc_matrix_B,
                                         const size_t dimension,
                                         Evaluator &evaluator,
                                         CKKSEncoder &encoder,
                                         GaloisKeys &galois_keys,
                                         RelinKeys &relin_keys, double scale) {
  size_t n_cols_B = enc_matrix_B.size();
  vector<Ciphertext> result(n_cols_B);

  for (size_t j = 0; j < n_cols_B; ++j) {
    result[j] = rowwise_encmat_encvec(enc_matrix_A, enc_matrix_B[j], dimension,
                                      evaluator, encoder, galois_keys,
                                      relin_keys, scale);
  }
  return result;
}

Ciphertext diagonal_encmat_encvec(const vector<Ciphertext> &enc_matrix_A,
                                  const Ciphertext &enc_vector_x,
                                  const size_t dimension, Evaluator &evaluator,
                                  CKKSEncoder &encoder, GaloisKeys &galois_keys,
                                  RelinKeys &relin_keys, double scale) {
  Ciphertext result;
  bool first = true;

  for (size_t d = 0; d < dimension; ++d) {
    Ciphertext rotated_vec;
    if (d > 0) {
      PROFILE_START("Rotate Vector");
      evaluator.rotate_vector(enc_vector_x, static_cast<int>(d), galois_keys,
                              rotated_vec);
      PROFILE_END("Rotate Vector");
    } else {
      rotated_vec = enc_vector_x;
    }

    Ciphertext prod;
    PROFILE_START("Homomorphic Mult");
    evaluator.multiply(enc_matrix_A[d], rotated_vec, prod);
    evaluator.relinearize_inplace(prod, relin_keys);
    evaluator.rescale_to_next_inplace(prod);
    PROFILE_END("Homomorphic Mult");

    if (first) {
      result = prod;
      first = false;
    } else {
      PROFILE_START("Result Add");
      evaluator.add_inplace(result, prod);
      PROFILE_END("Result Add");
    }
  }
  return result;
}

vector<Ciphertext>
diagonal_encmat_encmat(const vector<Ciphertext> &enc_matrix_A,
                       const vector<Ciphertext> &enc_matrix_B,
                       const size_t dimension, Evaluator &evaluator,
                       CKKSEncoder &encoder, GaloisKeys &galois_keys,
                       RelinKeys &relin_keys, double scale) {
  size_t n_cols_B = enc_matrix_B.size();
  vector<Ciphertext> result(n_cols_B);

  for (size_t j = 0; j < n_cols_B; ++j) {
    result[j] = diagonal_encmat_encvec(enc_matrix_A, enc_matrix_B[j], dimension,
                                       evaluator, encoder, galois_keys,
                                       relin_keys, scale);
  }
  return result;
}
} // namespace matvec
#endif // MATVEC_HPP
