// Sequential code: Qubits vs Execution Time
#include "quest.h"
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <iostream>

using qreal = qreal;

const int NUM_LAYERS = 20;  // Fixed number of layers
const qreal FIXED_NOISE = 0.05;  // Fixed noise level

// Test different qubit counts
std::vector<int> QUBIT_COUNTS = {4, 6, 8, 10, 12, 14, 16, 18, 20};

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

int main() {
    srand(time(0));
    initQuESTEnv();

    std::ofstream outFile("qubits_vs_time_sequential.csv");
    outFile << "NumQubits,ExecutionTime_sec\n";

    std::cout << "Running Sequential Qubits vs Time Benchmark...\n";
    std::cout << "Fixed parameters: " << NUM_LAYERS << " layers, " 
              << FIXED_NOISE << " noise level\n\n";

    for (int numQubits : QUBIT_COUNTS) {
        std::cout << "Testing with " << numQubits << " qubits... " << std::flush;

        int totalParams = numQubits * 2 * NUM_LAYERS;
        std::vector<qreal> angles = generateRandomAngles(totalParams);

        Qureg qureg = createQureg(numQubits);
        initZeroState(qureg);

        // Measure execution time
        double startTime = omp_get_wtime();

        for (int layer = 0; layer < NUM_LAYERS; ++layer) {
            applyHEALayer(qureg, angles, layer, numQubits);
        }

        double endTime = omp_get_wtime();
        double executionTime = endTime - startTime;

        std::cout << executionTime << " seconds\n";

        outFile << numQubits << "," << executionTime << "\n";
        outFile.flush();

        destroyQureg(qureg);
    }

    outFile.close();

    std::cout << "\nBenchmark complete! Results saved to qubits_vs_time_sequential.csv\n";
    return 0;
}
