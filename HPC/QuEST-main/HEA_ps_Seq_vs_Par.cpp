// HEA_ps_Seq_vs_Par.cpp
#include "quest.h"

#include <cstdio>     // For printf
#include <fstream>    // For ofstream
#include <vector>     // For std::vector
#include <omp.h>      // For OpenMP functions like omp_get_wtime
#include <cstdlib>    // For srand, rand
#include <ctime>      // For time

using qreal = qreal;

const qreal FIXED_NOISE = 0.5; // Small parameter noise level

// Function Declarations
std::vector<qreal> generateRandomAngles(int numParams);
void applyNoisyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex, qreal noise);

int main() {
    srand(time(0));
    initQuESTEnv();

    // Define ranges
    const std::vector<int> QUBIT_RANGE = {10, 15, 20, 25, 30};
    const std::vector<int> LAYER_RANGE = {10, 25, 50, 75, 100};

    const int NUM_RUNS = 5; // Run multiple simulations per thread

    std::ofstream resultFile("benchmark_results.csv");
    if (!resultFile.is_open()) {
        printf("Error: Could not open benchmark_results.csv\n");
        return -1;
    }

    resultFile << "Qubits,Layers,Mode,Time,AvgError\n";

    for (int NUM_QUBITS : QUBIT_RANGE) {
        for (int NUM_LAYERS : LAYER_RANGE) {

            int totalParams = NUM_QUBITS * 2 * NUM_LAYERS;
            std::vector<qreal> fixedAngles = generateRandomAngles(totalParams);

            // Sequential Execution
            qreal seq_total_error = 0.0;
            double seq_start = omp_get_wtime();

            for (int run = 0; run < NUM_RUNS; ++run) {
                Qureg qureg_seq = createQureg(NUM_QUBITS);
                initZeroState(qureg_seq);

                for (int layer = 0; layer < NUM_LAYERS; ++layer) {
                    applyNoisyHEALayer(qureg_seq, fixedAngles, layer, FIXED_NOISE);
                }

                Qureg target_seq = createQureg(NUM_QUBITS);
                initZeroState(target_seq);

                qreal fidelity = calcFidelity(qureg_seq, target_seq);
                seq_total_error += (1.0 - fidelity);

                destroyQureg(qureg_seq);
                destroyQureg(target_seq);
            }

            double seq_end = omp_get_wtime();
            double seq_time = seq_end - seq_start;
            qreal seq_avg_error = seq_total_error / NUM_RUNS;

            // Parallel Execution
            qreal par_total_error = 0.0;
            double par_start = omp_get_wtime();

            #pragma omp parallel reduction(+:par_total_error)
            {
                #pragma omp for schedule(static)
                for (int run = 0; run < NUM_RUNS; ++run) {
                    Qureg qureg_par = createQureg(NUM_QUBITS);
                    initZeroState(qureg_par);

                    for (int layer = 0; layer < NUM_LAYERS; ++layer) {
                        applyNoisyHEALayer(qureg_par, fixedAngles, layer, FIXED_NOISE);
                    }

                    Qureg target_par = createQureg(NUM_QUBITS);
                    initZeroState(target_par);

                    qreal fidelity = calcFidelity(qureg_par, target_par);
                    par_total_error += (1.0 - fidelity);

                    destroyQureg(qureg_par);
                    destroyQureg(target_par);
                }
            }

            double par_end = omp_get_wtime();
            double par_time = par_end - par_start;
            qreal par_avg_error = par_total_error / NUM_RUNS;

            // Write to CSV
            resultFile << NUM_QUBITS << "," << NUM_LAYERS << ",Sequential," << seq_time << "," << seq_avg_error << "\n";
            resultFile << NUM_QUBITS << "," << NUM_LAYERS << ",Parallel," << par_time << "," << par_avg_error << "\n";

            // Print progress
            printf("Q=%2d L=%3d | Seq=%.4fs Par=%.4fs Speedup=%.2fx Error=%.6f\n",
                   NUM_QUBITS, NUM_LAYERS,
                   seq_time, par_time, seq_time / par_time, par_avg_error);
        }
    }

    resultFile.close();
    printf("Benchmark complete. Results saved to benchmark_results.csv\n");

    return 0;
}

// Generate random angles between 0 and pi
std::vector<qreal> generateRandomAngles(int numParams) {
    std::vector<qreal> angles(numParams);
    for (int i = 0; i < numParams; ++i)
        angles[i] = ((qreal)rand() / RAND_MAX) * M_PI;
    return angles;
}

// Apply one noisy HEA layer
void applyNoisyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex, qreal noise) {
    int numQubits = qureg.numQubits;
    int offset = layerIndex * numQubits * 2;

    for (int qubit = 0; qubit < numQubits; ++qubit) {
        qreal rx_angle = angles[offset + qubit * 2];
        qreal ry_angle = angles[offset + qubit * 2 + 1];

        rx_angle *= (1 + noise * ((qreal)rand() / RAND_MAX - 0.5));
        ry_angle *= (1 + noise * ((qreal)rand() / RAND_MAX - 0.5));

        applyRotateX(qureg, qubit, rx_angle);
        applyRotateY(qureg, qubit, ry_angle);
    }

    for (int qubit = 0; qubit < numQubits - 1; ++qubit) {
        applyControlledMultiQubitNot(qureg, qubit, new int[1]{qubit + 1}, 1);
    }
}
