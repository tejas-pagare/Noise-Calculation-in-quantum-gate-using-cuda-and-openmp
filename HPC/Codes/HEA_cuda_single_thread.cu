// CUDA Single-Thread HEA Simulation with Qubits vs Time Analysis
// This code runs the HEA simulation on GPU using QuEST with CUDA backend
// and generates data for plotting Qubits vs Execution Time graph
#include "quest.h"
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <chrono>
#include <iomanip>

using qreal = qreal;

// Configuration - adjust as needed
const int NUM_LAYERS = 10;          // Number of HEA layers
const qreal NOISE_LEVEL = 0.2;      // Noise probability

// Range of qubits to test (for density matrix, keep max reasonable for GPU memory)
const int MIN_QUBITS = 2;
const int MAX_QUBITS = 14;          // Adjust based on GPU memory
const int QUBIT_STEP = 1;

std::vector<qreal> generateRandomAngles(int numParams) {
    std::vector<qreal> angles(numParams);
    for (int i = 0; i < numParams; ++i)
        angles[i] = ((qreal)rand() / RAND_MAX) * M_PI;
    return angles;
}

void applyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex, int numQubits) {
    int offset = layerIndex * numQubits * 2;

    // Single-qubit rotations
    for (int qubit = 0; qubit < numQubits; ++qubit) {
        applyRotateX(qureg, qubit, angles[offset + qubit * 2]);
        applyRotateY(qureg, qubit, angles[offset + qubit * 2 + 1]);
    }

    // Entanglement layer (linear chain of CNOTs)
    for (int qubit = 0; qubit < numQubits - 1; ++qubit) {
        int targets[] = {qubit + 1};
        applyControlledMultiQubitNot(qureg, qubit, targets, 1);
    }
}

void applyDepolarizingNoise(Qureg qureg, int qubit, qreal probability) {
    // Apply depolarizing noise using damping channel
    mixDamping(qureg, qubit, probability);
}

void simulateNoisyHEA(Qureg qureg, const std::vector<qreal>& angles, qreal noiseProb, int numQubits, int numLayers) {
    for (int layer = 0; layer < numLayers; ++layer) {
        // Apply HEA layer
        applyHEALayer(qureg, angles, layer, numQubits);
        
        // Apply noise to all qubits
        for (int qubit = 0; qubit < numQubits; ++qubit) {
            applyDepolarizingNoise(qureg, qubit, noiseProb);
        }
    }
}

double runSimulationForQubits(int numQubits, int numLayers, qreal noiseLevel) {
    // Generate random angles for HEA circuit
    int totalParams = numQubits * 2 * numLayers;
    std::vector<qreal> angles = generateRandomAngles(totalParams);

    // Create quantum register as DENSITY MATRIX for noise simulation
    Qureg qureg = createDensityQureg(numQubits);
    initZeroState(qureg);

    // Time the simulation
    auto start = std::chrono::high_resolution_clock::now();
    
    simulateNoisyHEA(qureg, angles, noiseLevel, numQubits, numLayers);
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // Cleanup
    destroyQureg(qureg);

    return elapsed.count();
}

int main() {
    srand(time(0));
    
    // Initialize QuEST environment (will use GPU if compiled with CUDA)
    initQuESTEnv();

    std::cout << "============================================================\n";
    std::cout << "  CUDA Single-Thread HEA Simulation - Qubits vs Time Analysis\n";
    std::cout << "============================================================\n";
    std::cout << "Configuration:\n";
    std::cout << "  Qubit Range: " << MIN_QUBITS << " to " << MAX_QUBITS << "\n";
    std::cout << "  Layers: " << NUM_LAYERS << "\n";
    std::cout << "  Noise Level: " << NOISE_LEVEL << "\n";
    std::cout << "============================================================\n\n";

    // Open CSV file for results
    std::ofstream csvFile("cuda_qubits_vs_time.csv");
    csvFile << "Qubits,Time(seconds),DensityMatrixSize,GPUMemoryEstimate(MB)\n";

    // Store results for display
    std::vector<int> qubitCounts;
    std::vector<double> executionTimes;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n" << std::setw(10) << "Qubits" 
              << std::setw(18) << "Time (sec)" 
              << std::setw(22) << "Matrix Size"
              << std::setw(20) << "Est. Memory (MB)" << "\n";
    std::cout << std::string(70, '-') << "\n";

    for (int numQubits = MIN_QUBITS; numQubits <= MAX_QUBITS; numQubits += QUBIT_STEP) {
        // Calculate density matrix size: 2^n x 2^n complex numbers
        long long matrixDim = 1LL << numQubits;  // 2^n
        long long totalElements = matrixDim * matrixDim;
        double memoryMB = (totalElements * 2 * sizeof(qreal)) / (1024.0 * 1024.0);  // 2 for complex

        std::cout << std::setw(10) << numQubits << std::flush;

        try {
            double time = runSimulationForQubits(numQubits, NUM_LAYERS, NOISE_LEVEL);
            
            qubitCounts.push_back(numQubits);
            executionTimes.push_back(time);

            std::cout << std::setw(18) << time 
                      << std::setw(22) << (std::to_string(matrixDim) + "x" + std::to_string(matrixDim))
                      << std::setw(20) << memoryMB << "\n";

            // Write to CSV
            csvFile << numQubits << "," << time << "," << totalElements << "," << memoryMB << "\n";
            csvFile.flush();
        }
        catch (...) {
            std::cout << std::setw(18) << "FAILED (OOM?)" 
                      << std::setw(22) << (std::to_string(matrixDim) + "x" + std::to_string(matrixDim))
                      << std::setw(20) << memoryMB << "\n";
            break;  // Stop if we run out of memory
        }
    }

    csvFile.close();

    // Print summary
    std::cout << "\n============================================================\n";
    std::cout << "  Results Summary\n";
    std::cout << "============================================================\n";
    std::cout << "Data saved to: cuda_qubits_vs_time.csv\n";
    std::cout << "Total configurations tested: " << qubitCounts.size() << "\n";
    
    if (!executionTimes.empty()) {
        double minTime = *std::min_element(executionTimes.begin(), executionTimes.end());
        double maxTime = *std::max_element(executionTimes.begin(), executionTimes.end());
        std::cout << "Time range: " << minTime << " - " << maxTime << " seconds\n";
    }
    std::cout << "============================================================\n";

    // Generate Python plotting script
    std::ofstream plotScript("plot_cuda_qubits_vs_time.py");
    plotScript << R"(#!/usr/bin/env python3
"""
Plot Qubits vs Execution Time for CUDA HEA Simulation
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read data
df = pd.read_csv('cuda_qubits_vs_time.csv')

# Create figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Linear scale
ax1 = axes[0]
ax1.plot(df['Qubits'], df['Time(seconds)'], 'b-o', linewidth=2, markersize=8, label='CUDA (Single Thread)')
ax1.set_xlabel('Number of Qubits', fontsize=12)
ax1.set_ylabel('Execution Time (seconds)', fontsize=12)
ax1.set_title('Qubits vs Execution Time (Linear Scale)', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.legend()

# Plot 2: Log scale (shows exponential scaling)
ax2 = axes[1]
ax2.semilogy(df['Qubits'], df['Time(seconds)'], 'r-s', linewidth=2, markersize=8, label='CUDA (Single Thread)')
ax2.set_xlabel('Number of Qubits', fontsize=12)
ax2.set_ylabel('Execution Time (seconds) - Log Scale', fontsize=12)
ax2.set_title('Qubits vs Execution Time (Log Scale)', fontsize=14)
ax2.grid(True, alpha=0.3, which='both')
ax2.legend()

plt.tight_layout()
plt.savefig('cuda_qubits_vs_time.png', dpi=150, bbox_inches='tight')
plt.savefig('cuda_qubits_vs_time.pdf', bbox_inches='tight')
print("Plots saved as 'cuda_qubits_vs_time.png' and 'cuda_qubits_vs_time.pdf'")

# Show plot
plt.show()

# Print statistics
print("\n=== Statistics ===")
print(f"Qubits tested: {df['Qubits'].min()} to {df['Qubits'].max()}")
print(f"Min time: {df['Time(seconds)'].min():.6f} sec at {df.loc[df['Time(seconds)'].idxmin(), 'Qubits']} qubits")
print(f"Max time: {df['Time(seconds)'].max():.6f} sec at {df.loc[df['Time(seconds)'].idxmax(), 'Qubits']} qubits")

# Calculate scaling factor (should be ~4x per qubit for density matrix)
if len(df) > 1:
    scaling_factors = []
    for i in range(1, len(df)):
        if df['Time(seconds)'].iloc[i-1] > 0:
            factor = df['Time(seconds)'].iloc[i] / df['Time(seconds)'].iloc[i-1]
            scaling_factors.append(factor)
    if scaling_factors:
        avg_scaling = np.mean(scaling_factors)
        print(f"Average scaling factor per qubit: {avg_scaling:.2f}x")
        print("(Expected ~4x for density matrix simulation)")
)";
    plotScript.close();

    std::cout << "\nPython plotting script saved as: plot_cuda_qubits_vs_time.py\n";
    std::cout << "Run with: python3 plot_cuda_qubits_vs_time.py\n";

    return 0;
}
