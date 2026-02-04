#include "quest.h"
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>

using qreal = qreal;

const int NUM_QUBITS = 10;
const int NUM_LAYERS = 50;

std::vector<qreal> generateRandomAngles(int numParams) {
    std::vector<qreal> angles(numParams);
    for (int i = 0; i < numParams; ++i)
        angles[i] = ((qreal)rand() / RAND_MAX) * M_PI;
    return angles;
}

void applyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex, qreal noise) {
    int offset = layerIndex * NUM_QUBITS * 2;

    for (int qubit = 0; qubit < NUM_QUBITS; ++qubit) {
        applyRotateX(qureg, qubit, angles[offset + qubit * 2]);
        applyRotateY(qureg, qubit, angles[offset + qubit * 2 + 1]);
        qreal error;
        mixDephasing(qureg, qubit, noise, &error);
    }

    for (int qubit = 0; qubit < NUM_QUBITS - 1; ++qubit) {
        applyControlledMultiQubitNot(qureg, qubit, new int[1]{qubit + 1}, 1);
    }
}

int main() {
    srand(time(0));
    initQuESTEnv();

    int totalParams = NUM_QUBITS * 2 * NUM_LAYERS;
    std::vector<qreal> angles = generateRandomAngles(totalParams);

    Qureg qureg_seq = createDensityQureg(NUM_QUBITS);
    double seq_time = 0.0;

    {
        double start = omp_get_wtime();
        for (int rep = 0; rep < 5; ++rep) {
            initZeroState(qureg_seq);
            for (int layer = 0; layer < NUM_LAYERS; ++layer) {
                applyHEALayer(qureg_seq, angles, layer, 0.10);
            }
        }
        double end = omp_get_wtime();
        seq_time = (end - start) / 5.0;
    }

    destroyQureg(qureg_seq);

    double par_time = 0.0;
    double par_start = omp_get_wtime();

    #pragma omp parallel for
    for (int rep = 0; rep < 5; ++rep) {
        Qureg qureg_par = createDensityQureg(NUM_QUBITS);
        initZeroState(qureg_par);
        for (int layer = 0; layer < NUM_LAYERS; ++layer) {
            applyHEALayer(qureg_par, angles, layer, 0.10);
        }
        destroyQureg(qureg_par);
    }

    par_time = (omp_get_wtime() - par_start) / 5.0;

    std::ofstream resultFile("seq_vs_parallel.csv");
    resultFile << "Mode,Time,Speedup\n";
    resultFile << "Sequential," << seq_time << ",1.00\n";
    resultFile << "Parallel," << par_time << "," << (seq_time / par_time) << "\n";

    resultFile.close();
    printf("Sequential Time: %.4f s\n", seq_time);
    printf("Parallel Time:   %.4f s\n", par_time);
    printf("Speedup: %.2fx\n", seq_time / par_start);

    return 0;
}
