#include "quest.h"
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>

using qreal = qreal;

// Configuration
const int NUM_QUBITS = 23;         // Simulate 10-qubit HEA
const int NUM_LAYERS = 50;        // Number of HEA layers

std::vector<qreal> NOISE_LEVELS = {
    0.00, 0.02, 0.04, 0.06,
    0.08, 0.10, 0.12, 0.14,
    0.16, 0.18, 0.20, 0.22,
    0.24, 0.26, 0.28, 0.30,
    0.32, 0.34, 0.36, 0.38,
    0.40, 0.42, 0.44, 0.46,
    0.48, 0.50
};

// Generate random angles for HEA
std::vector<qreal> generateRandomAngles(int numParams) {
    std::vector<qreal> angles(numParams);
    for (int i = 0; i < numParams; ++i)
        angles[i] = ((qreal)rand() / RAND_MAX) * M_PI;
    return angles;
}

// Apply one layer of HEA on N qubits
void applyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex) {
    int offset = layerIndex * NUM_QUBITS * 2;

    // Single-qubit rotations
    for (int qubit = 0; qubit < NUM_QUBITS; ++qubit) {
        applyRotateX(qureg, qubit, angles[offset + qubit * 2]);       // RX
        applyRotateY(qureg, qubit, angles[offset + qubit * 2 + 1]);   // RY
    }

    // Entanglement layer (linear chain of CNOTs)
    for (int qubit = 0; qubit < NUM_QUBITS - 1; ++qubit) {
        applyControlledMultiQubitNot(qureg, qubit, new int[1]{qubit + 1}, 1);
    }
}

// Simulate HEA and measure execution time
double simulateHEA(Qureg qureg, const std::vector<qreal>& angles) {
    initZeroState(qureg);

    double start = omp_get_wtime();

    for (int layer = 0; layer < NUM_LAYERS; ++layer) {
        applyHEALayer(qureg, angles, layer);
    }

    double end = omp_get_wtime();
    return end - start;
}

int main() {
    srand(time(0)); // Seed randomness

    initQuESTEnv();
    QuESTEnv env = getQuESTEnv();

    int totalParams = NUM_QUBITS * 2 * NUM_LAYERS;
    std::vector<qreal> fixedAngles = generateRandomAngles(totalParams);

    // Run sequential version
    Qureg qureg_seq = createQureg(NUM_QUBITS);
    double seq_total_time = 0.0;

    {
        double start = omp_get_wtime();
        for (size_t i = 0; i < NOISE_LEVELS.size(); ++i) {
            initZeroState(qureg_seq);
            for (int layer = 0; layer < NUM_LAYERS; ++layer) {
                applyHEALayer(qureg_seq, fixedAngles, layer);
            }
        }
        double end = omp_get_wtime();
        seq_total_time = end - start;
    }

    destroyQureg(qureg_seq);

    // Run parallel version
    double par_total_time = 0.0;
    double par_start = omp_get_wtime();

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < NOISE_LEVELS.size(); ++i) {
        qreal noise = NOISE_LEVELS[i];
        Qureg qureg_local = createQureg(NUM_QUBITS);
        initZeroState(qureg_local);

        for (int layer = 0; layer < NUM_LAYERS; ++layer) {
            applyHEALayer(qureg_local, fixedAngles, layer);
        }

        destroyQureg(qureg_local);
    }

    par_total_time = omp_get_wtime() - par_start;

    // Output summary to CSV (with metric)
    std::ofstream resultFile("timing_summary.csv", std::ios::app);

    if (resultFile.tellp() == 0) {
        // Write header only once
        resultFile << "NoiseLevel,SequentialTime,ParallelTime\n";
    }

    // Write per-noise-level lines (optional)
    for (auto noise : NOISE_LEVELS) {
        resultFile << noise << "," << seq_total_time / NOISE_LEVELS.size() << "," << par_total_time / NOISE_LEVELS.size() << "\n";
    }

    // Write final metric line
    resultFile << "#TotalSequential,TotalParallel,Speedup\n";
    resultFile << "#" << seq_total_time << "," << par_total_time << "," << (seq_total_time / par_total_time) << "\n";

    resultFile.close();

    // Also print to console
    printf("Total Sequential Time: %.6f seconds\n", seq_total_time);
    printf("Total Parallel Time: %.6f seconds\n", par_total_time);
    printf("Speedup: %.2fx\n", seq_total_time / par_total_time);

    return 0;
}
