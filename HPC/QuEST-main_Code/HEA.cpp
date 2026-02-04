#include "quest.h"
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>

using qreal = qreal;

// Configuration
const int NUM_QUBITS = 8;         // Simulate 4-qubit HEA
const int NUM_LAYERS = 25;         // Number of HEA layers
std::vector<qreal> NOISE_LEVELS = {
    0.00, 0.02, 0.04, 0.06,
    0.08, 0.10, 0.12, 0.14,
    0.16, 0.18, 0.20, 0.22,
    0.24, 0.26, 0.28, 0.30
};

// Generate random angles for HEA
std::vector<qreal> generateRandomAngles(int numParams) {
    std::vector<qreal> angles(numParams);
    for (int i = 0; i < numParams; ++i)
        angles[i] = ((qreal)rand() / RAND_MAX) * M_PI;
    return angles;
}

// Apply one layer of HEA on N qubits
void applyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex,
                   std::ofstream& traceFile, double startTime) {
    int offset = layerIndex * NUM_QUBITS * 2;

    // Single-qubit rotations
    for (int qubit = 0; qubit < NUM_QUBITS; ++qubit) {
        applyRotateX(qureg, qubit, angles[offset + qubit * 2]);       // RX
	applyRotateY(qureg, qubit, angles[offset + qubit * 2 + 1]);   // RY

        double t = omp_get_wtime();
        traceFile << t - startTime << "," << qubit << ",RX," << angles[offset + qubit * 2] << "\n";
        traceFile << t - startTime << "," << qubit << ",RY," << angles[offset + qubit * 2 + 1] << "\n";
    }

    // Entanglement layer (linear chain of CNOTs)
    for (int qubit = 0; qubit < NUM_QUBITS - 1; ++qubit) {
        applyControlledMultiQubitNot(qureg, qubit, new int[1]{qubit + 1}, 1);

        double t = omp_get_wtime();
        traceFile << t - startTime << "," << qubit << ",CNOT," << qubit + 1 << "\n";
    }
}

// Simulate HEA with noise and track time
void simulateHEA(Qureg qureg, qreal noise, bool useTracked, std::ofstream& outFile, std::ofstream& traceFile, double& totalTime) {
    initZeroState(qureg);

    int totalParams = NUM_QUBITS * 2 * NUM_LAYERS;
    std::vector<qreal> angles = generateRandomAngles(totalParams);

    double start = omp_get_wtime();

    std::vector<double> layerTimes(NUM_LAYERS, 0.0);

    // Apply multiple layers with timing per layer
    for (int layer = 0; layer < NUM_LAYERS; ++layer) {
        double layerStart = omp_get_wtime();

        applyHEALayer(qureg, angles, layer, traceFile, start);

        // Apply dephasing noise to all qubits after each layer
        qreal totalError = 0;
        for (int qubit = 0; qubit < NUM_QUBITS; ++qubit) {
            qreal error = 0;
            if (useTracked)
                mixDephasingTracked(qureg, qubit, noise, &error);
            else
                mixDephasing(qureg, qubit, noise);
            totalError += error;
        }

        double layerEnd = omp_get_wtime();
        layerTimes[layer] = layerEnd - layerStart;
    }

    double end = omp_get_wtime();
    totalTime = end - start;

    qreal totalError = 0;
for (int qubit = 0; qubit < NUM_QUBITS; ++qubit) {
    qreal error = 0;
    if (useTracked)
        mixDephasingTracked(qureg, qubit, noise, &error);
    else
        mixDephasing(qureg, qubit, noise);
    totalError += error;
}
qreal avgError = totalError / NUM_QUBITS;

    // Output results
    outFile << "HEA," << NUM_QUBITS << "," << NUM_LAYERS << "," << noise << ","
            << useTracked << "," << avgError << "\n";

    // Print layer-wise times
    printf("Noise=%.2f, Tracking=%d\n", noise, useTracked);
    for (int layer = 0; layer < NUM_LAYERS; ++layer) {
        printf("  Layer %2d: %.6f s\n", layer, layerTimes[layer]);
    }
    printf("  Total: %.6f s\n\n", totalTime);
}

int main() {
    srand(time(0)); // Seed randomness

    initQuESTEnv();
    QuESTEnv env = getQuESTEnv();

    Qureg qureg = createDensityQureg(NUM_QUBITS);

    // Open output files
    std::ofstream out("hea_results.csv");
    out << "Circuit,NumQubits,NumLayers,NoiseProbability,UseTracked,MeasuredError\n";

    std::ofstream traceFile("hea_trace.csv");
    traceFile << "Time,Qubit,Gate,TargetQubit,Angle\n";

    double totalStart = omp_get_wtime();

    // Run simulation for all noise levels and tracking modes
    for (auto noise : NOISE_LEVELS) {
        double timeTrack = 0, timeNoTrack = 0;

        simulateHEA(qureg, noise, true, out, traceFile, timeTrack);   // With tracking
        simulateHEA(qureg, noise, false, out, traceFile, timeNoTrack); // Without tracking

        printf("Noise=%.2f => Track: %.6f s, NoTrack: %.6f s\n\n", noise, timeTrack, timeNoTrack);
    }

    double totalEnd = omp_get_wtime();
    printf("Total Simulation Time: %.6f seconds\n", totalEnd - totalStart);

    destroyQureg(qureg);

    out.close();
    traceFile.close();

    return 0;
}
