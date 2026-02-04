// qubits_layers_benchmark.cpp
#include "quest.h"
#include <cstdio>
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>

using qreal = qreal;

const qreal FIXED_NOISE = 0.5; // Noise level applied to RX/RY angles

// Define sweep ranges
const std::vector<int> QUBIT_RANGE = {4, 6, 8, 10};     // Safe for density matrices
const std::vector<int> LAYER_RANGE = {10, 25, 50};       // Test different depths

std::vector<qreal> generateRandomAngles(int numParams) {
    std::vector<qreal> angles(numParams);
    for (int i = 0; i < numParams; ++i)
        angles[i] = ((qreal)rand() / RAND_MAX) * M_PI;
    return angles;
}

void applyNoisyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex, qreal noise) {
    int numQubits = qureg.numQubits;
    int offset = layerIndex * numQubits * 2;

    for (int qubit = 0; qubit < numQubits; ++qubit) {
        qreal rx_angle = angles[offset + qubit * 2];
        qreal ry_angle = angles[offset + qubit * 2 + 1];

        if (noise > 0) {
            rx_angle *= (1 + noise * ((qreal)rand() / RAND_MAX - 0.5));
            ry_angle *= (1 + noise * ((qreal)rand() / RAND_MAX - 0.5));
        }

        applyRotateX(qureg, qubit, rx_angle);
        applyRotateY(qureg, qubit, ry_angle);
    }

    for (int qubit = 0; qubit < numQubits - 1; ++qubit) {
        applyControlledMultiQubitNot(qureg, qubit, new int[1]{qubit + 1}, 1);
    }
}

int main() {
    srand(time(0));
    initQuESTEnv();

    // Set number of threads (based on threads_vs_time.csv: best at ~6 threads)
    omp_set_num_threads(6);

    std::ofstream resultFile("benchmark_optimized.csv");
    if (!resultFile.is_open()) {
        printf("Error: Could not open output file\n");
        return -1;
    }

    resultFile << "Qubits,Layers,SeqTime(s),ParTime(s),Speedup,AvgError\n";

    // Flatten qubit-layer combinations for full OpenMP scheduling
    std::vector<std::pair<int, int>> combos;
    for (int q : QUBIT_RANGE)
        for (int l : LAYER_RANGE)
            combos.push_back({q, l});

    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (size_t combo_idx = 0; combo_idx < combos.size(); ++combo_idx) {
            int NUM_QUBITS = combos[combo_idx].first;
            int NUM_LAYERS = combos[combo_idx].second;
            int totalParams = NUM_QUBITS * 2 * NUM_LAYERS;

            std::vector<qreal> angles = generateRandomAngles(totalParams);

            // Sequential Execution
            double seq_start = omp_get_wtime();
            {
                Qureg noisyReg_seq = createQureg(NUM_QUBITS);
                Qureg cleanReg_seq = createQureg(NUM_QUBITS);

                initZeroState(noisyReg_seq);
                initZeroState(cleanReg_seq);

                for (int layer = 0; layer < NUM_LAYERS; ++layer) {
                    applyNoisyHEALayer(noisyReg_seq, angles, layer, FIXED_NOISE);
                    applyNoisyHEALayer(cleanReg_seq, angles, layer, 0.0); // No noise
                }

                qreal fidelity = calcFidelity(noisyReg_seq, cleanReg_seq);
                qreal error = 1.0 - fidelity;

                destroyQureg(noisyReg_seq);
                destroyQureg(cleanReg_seq);

                double seq_end = omp_get_wtime();
                double seq_time = seq_end - seq_start;

                // Parallel Execution
                double par_start = omp_get_wtime();
                {
                    Qureg noisyReg_par = createQureg(NUM_QUBITS);
                    Qureg cleanReg_par = createQureg(NUM_QUBITS);

                    initZeroState(noisyReg_par);
                    initZeroState(cleanReg_par);

                    for (int layer = 0; layer < NUM_LAYERS; ++layer) {
                        applyNoisyHEALayer(noisyReg_par, angles, layer, FIXED_NOISE);
                        applyNoisyHEALayer(cleanReg_par, angles, layer, 0.0); // No noise
                    }

                    qreal fidelity = calcFidelity(noisyReg_par, cleanReg_par);
                    qreal par_error = 1.0 - fidelity;

                    destroyQureg(noisyReg_par);
                    destroyQureg(cleanReg_par);

                    double par_end = omp_get_wtime();
                    double par_time = par_end - par_start;

                    // Write to CSV safely
                    #pragma omp critical
                    {
                        qreal speedup = seq_time / par_time;

                        resultFile << NUM_QUBITS << "," 
                                   << NUM_LAYERS << ","
                                   << seq_time << "," 
                                   << par_time << ","
                                   << speedup << "," 
                                   << par_error << "\n";

                        printf("Q=%2d L=%3d | Seq=%.4fs Par=%.4fs Speedup=%.2fx Error=%.6f\n",
                               NUM_QUBITS, NUM_LAYERS,
                               seq_time, par_time, speedup, par_error);
                    }
                }
            }
        }
    }

    resultFile.close();
    printf("Benchmark complete. Results saved to benchmark_optimized.csv\n");

    return 0;
}
