#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace utils
{
    using namespace std;

    // Generating a matrix of a desired size by random values
    vector<vector<double>> random_matrix(
        size_t rows,
        size_t cols,
        double min_val = -10.0,
        double max_val = 10.0)
    {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(min_val, max_val);

        vector<vector<double>> mat(rows, vector<double>(cols));
        for (auto &row : mat)
            for (auto &val : row)
                val = dis(gen);

        return mat;
    }

    // Generating a vector of a desired size by random values
    vector<double> random_vector(
        size_t rows,
        double min_val = -10.0,
        double max_val = 10.0)
    {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(min_val, max_val);

        vector<double> vec(rows, 0.0);
        for (auto &row : vec)
            row = dis(gen);
        return vec;
    }

    // printing matrix
    void print(const vector<vector<double>>& matrix)
    {
        for (const auto &row : matrix)
        {
            for (const auto &col : row)
                if (abs(col) > 1e-6)
                    cout << fixed << setprecision(5) << col << "\t";
                else
                    cout << "0.00000\t";
            cout << endl;
        }
    }

    // printing vector
    void print(const vector<double>& vec)
    {
        for (const auto &val : vec)
            if (abs(val) > 1e-6)
                cout << fixed << setprecision(5) << val << endl;
            else
                 cout << "0.00000" << endl;
    }

    // calculating the expected results for matrix-vector product
    vector<double> expected_result(
        const vector<vector<double>> &matrix,
        const vector<double> &vec)
    {
        size_t rows = matrix.size();
        size_t cols = vec.size();
        vector<double> result(rows, 0.0);

        for (size_t i = 0; i < rows; ++i)
            for (size_t j = 0; j < cols; ++j)
                result[i] += matrix[i][j] * vec[j];

        return result;
    }

    // calculating the expected results for matrix-matrix product
    vector<vector<double>> expected_result(
        const vector<vector<double>> &A,
        const vector<vector<double>> &B)
    {
        size_t n = A.size();
        size_t m = B[0].size();
        size_t k = B.size();
        vector<vector<double>> result(n, vector<double>(m, 0.0));

        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < m; j++)
            {
                for (size_t t = 0; t < k; t++)
                {
                    result[i][j] += A[i][t] * B[t][j];
                }
            }
        }
        return result;
    }

    /**
     * @brief Checks if two vectors are approximately equal within a given tolerance.
     */
    void check_approximation(const vector<double>& result, const vector<double>& expected, double tolerance) {
        bool match = true;
        size_t check_size = expected.size(); 

        if (result.size() < check_size) {
            cout << "  [ERROR] Decoded result vector is smaller than the expected vector!" << endl;
            match = false;
        }

        if (match) {
            for (size_t i = 0; i < check_size; ++i) {
                if (abs(result[i] - expected[i]) > tolerance) {
                    cout << "  [FAIL] Mismatch at index " << i << ":" << endl;
                    cout << "    > Expected: " << fixed << setprecision(5) << expected[i] << endl;
                    cout << "    > Actual:   " << fixed << setprecision(5) << result[i] << endl;
                    cout << "    > Delta:    " << abs(result[i] - expected[i]) << endl;
                    match = false;
                    break; 
                }
            }
        }

        if (match) {
            cout << "  [PASS] All " << check_size << " values match within a tolerance of " << tolerance << "." << endl;
        }
    }

    /**
     * @brief Checks if two matrices are approximately equal within a given tolerance.
     */
    void check_matrix_approximation(
        const vector<vector<double>>& result, 
        const vector<vector<double>>& expected, 
        double tolerance) 
    {
        if (result.size() != expected.size() || (result.empty() && !expected.empty())) {
            cout << "  [ERROR] Matrix row counts do not match!" << endl;
            return;
        }

        for (size_t i = 0; i < expected.size(); ++i) {
            if (result[i].size() != expected[i].size()) {
                cout << "  [ERROR] Matrix column counts do not match on row " << i << "!" << endl;
                return;
            }
            for (size_t j = 0; j < expected[i].size(); ++j) {
                if (abs(result[i][j] - expected[i][j]) > tolerance) {
                    cout << "  [FAIL] Mismatch at index (" << i << ", " << j << "):" << endl;
                    cout << "    > Expected: " << fixed << setprecision(5) << expected[i][j] << endl;
                    cout << "    > Actual:   " << fixed << setprecision(5) << result[i][j] << endl;
                    cout << "    > Delta:    " << abs(result[i][j] - expected[i][j]) << endl;
                    return;
                }
            }
        }
        cout << "  [PASS] All " << expected.size() << "x" << (expected.empty() ? 0 : expected[0].size()) 
             << " values match within a tolerance of " << tolerance << "." << endl;
    }

} // namespace utils

#endif // UTILS_HPP
