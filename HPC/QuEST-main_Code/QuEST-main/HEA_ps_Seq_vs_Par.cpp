// seq_vs_parallel_with_noise.cpp
#include "quest.h"
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>

using qreal = qreal;

const int NUM_QUBITS = 20;
const int NUM_LAYERS = 50;
const qreal FIXED_NOISE = 0.50; // Fixed noise level

std::vector<qreal> generateRandomAngles(int numParams) {
    std::vector<qreal> angles(numParams);
    for (int i = 0; i < numParams; ++i)
        angles[i] = ((qreal)rand() / RAND_MAX) * M_PI;
    return angles;
}

void applyNoisyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex, qreal noise) {
    int offset = layerIndex * NUM_QUBITS * 2;

    for (int qubit = 0; qubit < NUM_QUBITS; ++qubit) {
        applyRotateX(qureg, qubit, angles[offset + qubit * 2]);
        applyRotateY(qureg, qubit, angles[offset + qubit * 2 + 1]);

        mixDephasing(qureg, qubit, noise); // Apply dephasing noise
    }

    for (int qubit = 0; qubit < NUM_QUBITS - 1; ++qubit) {
        applyControlledMultiQubitNot(qureg, qubit, new int[1]{qubit + 1}, 1);
    }
}

int main() {
    srand(time(0));
    initQuESTEnv();

    int totalParams = NUM_QUBITS * 2 * NUM_LAYERS;
    std::vector<qreal> fixedAngles = generateRandomAngles(totalParams);

    // Sequential Execution
    Qureg qureg_seq = createDensityQureg(NUM_QUBITS);
    double seq_start = omp_get_wtime();

    initZeroState(qureg_seq);
    for (int layer = 0; layer < NUM_LAYERS; ++layer) {
        applyNoisyHEALayer(qureg_seq, fixedAngles, layer, FIXED_NOISE);
    }

    double seq_end = omp_get_wtime();
    double seq_time = seq_end - seq_start;

    // Compute error
    Qureg target = createDensityQureg(NUM_QUBITS);
    initZeroState(target);
    qreal fidelity = calcFidelity(qureg_seq, target);
    qreal error = 1.0 - fidelity;
    destroyQureg(target);

    destroyQureg(qureg_seq);

    // Parallel Execution
    double par_start = omp_get_wtime();
    qreal total_error = 0.0;

    #pragma omp parallel reduction(+:total_error)
    {
        Qureg qureg_par = createDensityQureg(NUM_QUBITS);
        initZeroState(qureg_par);

        for (int layer = 0; layer < NUM_LAYERS; ++layer) {
            applyNoisyHEALayer(qureg_par, fixedAngles, layer, FIXED_NOISE);
        }

        Qureg target_par = createDensityQureg(NUM_QUBITS);
        initZeroState(target_par);
        qreal fidelity = calcFidelity(qureg_par, target_par);
        total_error += (1.0 - fidelity);
        destroyQureg(qureg_par);
        destroyQureg(target_par);
    }

    double par_end = omp_get_wtime();
    double par_time = par_end - par_start;

    // Average error across all threads
    qreal avg_error = total_error / omp_get_max_threads();

    // Output to CSV
    std::ofstream resultFile("seq_vs_parallel_with_noise.csv");
    resultFile << "Mode,Time,Error\n";
    resultFile << "Sequential," << seq_time << "," << error << "\n";
    resultFile << "Parallel," << par_time << "," << avg_error << "\n";
    resultFile.close();

    // Print to console
    printf("Sequential Time: %.4f s, Error: %.6f\n", seq_time, error);
    printf("Parallel Time:   %.4f s, Avg Error: %.6f\n", par_time, avg_error);
    printf("Speedup: %.2fx\n", seq_time / par_time);

    return 0;
}
