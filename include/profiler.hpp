#ifndef PROFILING_HPP
#define PROFILING_HPP

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <chrono>
#include <string>
#include <iomanip>

using namespace std;
using namespace chrono;

// Maps for profiling data
static unordered_map<string, high_resolution_clock::time_point> starts;
static unordered_map<string, nanoseconds> durations;
static unordered_map<string, size_t> counts;

// Start the timer
#define PROFILE_START(key) \
    starts[key] = high_resolution_clock::now();

// End the timer
#define PROFILE_END(key) \
    do { \
        auto end = high_resolution_clock::now(); \
        auto elapsed = duration_cast<nanoseconds>(end - starts[key]); \
        durations[key] += elapsed; \
        ++counts[key]; \
    } while(0)

/**
 * @brief Clears all stored profiling data (durations and counts).
 * Call this between separate benchmark runs to prevent data contamination.
 */
inline void clear_profiling() {
    durations.clear();
    counts.clear();
    starts.clear();
}

// Print all results
inline void print_profiling(const string& output_file = "../benchmarks/results/noName.txt") {
    ofstream out(output_file);

    out << left << setw(30) << "Operation"
        << right << setw(15) << "TotalTime(ms)"
        << setw(10) << "Calls"
        << setw(15) << "AvgTime(us)" << '\n';

    for (auto& [key, total_ns] : durations) {
        size_t calls = counts[key];
        double total_ms = total_ns.count() / 1e6;
        double avg_us = total_ns.count() / 1e3 / calls;

        out << left << setw(30) << key
            << right << setw(15) << fixed << setprecision(3) << total_ms
            << setw(10) << calls
            << setw(15) << fixed << setprecision(3) << avg_us << '\n';
    }

    out.close();
}


#endif // PROFILING_HPP
