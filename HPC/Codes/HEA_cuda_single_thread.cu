// CUDA Single-Thread HEA Simulation
// This code runs the HEA simulation on GPU using QuEST with CUDA backend
#include "quest.h"
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <chrono>

using qreal = qreal;

// Configuration - adjust as needed
const int NUM_QUBITS = 12;          // Number of qubits
const int NUM_LAYERS = 30;          // Number of HEA layers
const qreal NOISE_LEVEL = 0.2;      // Noise probability

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

void simulateNoisyHEA(Qureg qureg, const std::vector<qreal>& angles, qreal noiseProb) {
    for (int layer = 0; layer < NUM_LAYERS; ++layer) {
        // Apply HEA layer
        applyHEALayer(qureg, angles, layer, NUM_QUBITS);
        
        // Apply noise to all qubits
        for (int qubit = 0; qubit < NUM_QUBITS; ++qubit) {
            applyDepolarizingNoise(qureg, qubit, noiseProb);
        }
    }
}

int main() {
    srand(time(0));
    
    // Initialize QuEST environment (will use GPU if compiled with CUDA)
    initQuESTEnv();

    std::cout << "========================================\n";
    std::cout << "  CUDA Single-Thread HEA Simulation\n";
    std::cout << "========================================\n";
    std::cout << "Configuration:\n";
    std::cout << "  Qubits: " << NUM_QUBITS << "\n";
    std::cout << "  Layers: " << NUM_LAYERS << "\n";
    std::cout << "  Noise Level: " << NOISE_LEVEL << "\n";
    std::cout << "========================================\n\n";

    // Generate random angles for HEA circuit
    int totalParams = NUM_QUBITS * 2 * NUM_LAYERS;
    std::vector<qreal> angles = generateRandomAngles(totalParams);

    // Create quantum register as DENSITY MATRIX for noise simulation
    // QuEST will automatically use GPU memory when compiled with CUDA
    Qureg qureg = createDensityQureg(NUM_QUBITS);
    initZeroState(qureg);

    std::cout << "Starting CUDA-accelerated simulation...\n";
    auto start = std::chrono::high_resolution_clock::now();

    // Run the noisy HEA simulation
    simulateNoisyHEA(qureg, angles, NOISE_LEVEL);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "\n========================================\n";
    std::cout << "  Results\n";
    std::cout << "========================================\n";
    std::cout << "Simulation Time: " << elapsed.count() << " seconds\n";
    std::cout << "========================================\n";

    // Cleanup
    destroyQureg(qureg);

    return 0;
}
