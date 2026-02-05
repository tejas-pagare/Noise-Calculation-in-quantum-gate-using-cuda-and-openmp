// Sequential HEA for gprof profiling
#include "quest.h"
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <chrono>

using qreal = qreal;

// Configuration for profiling
const int NUM_QUBITS = 12;          // Medium size for meaningful profiling
const int NUM_LAYERS = 30;          // Enough layers to show bottlenecks
const qreal NOISE_LEVEL = 0.2;

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
        applyControlledMultiQubitNot(qureg, qubit, new int[1]{qubit + 1}, 1);
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
    initQuESTEnv();

    std::cout << "Starting Sequential HEA Profiling Run...\n";
    std::cout << "Configuration:\n";
    std::cout << "  Qubits: " << NUM_QUBITS << "\n";
    std::cout << "  Layers: " << NUM_LAYERS << "\n";
    std::cout << "  Noise: " << NOISE_LEVEL << "\n\n";

    // Generate random angles
    int totalParams = NUM_QUBITS * 2 * NUM_LAYERS;
    std::vector<qreal> angles = generateRandomAngles(totalParams);

    // Create quantum register as DENSITY MATRIX for noise simulation
    Qureg qureg = createDensityQureg(NUM_QUBITS);
    initZeroState(qureg);

    auto start = std::chrono::high_resolution_clock::now();

    // Run the simulation
    simulateNoisyHEA(qureg, angles, NOISE_LEVEL);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Simulation completed in " << elapsed.count() << " seconds\n";
    std::cout << "gmon.out file will be generated for gprof analysis\n";

    destroyQureg(qureg);

    return 0;
}
