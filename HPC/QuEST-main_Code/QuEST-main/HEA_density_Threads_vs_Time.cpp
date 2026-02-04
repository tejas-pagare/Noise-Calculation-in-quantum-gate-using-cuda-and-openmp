#include "quest.h"
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>

using qreal = qreal;

const int NUM_QUBITS = 10;
const int NUM_LAYERS = 50;
const qreal FIXED_NOISE = 0.10;
const std::vector<int> THREAD_COUNTS = {1, 2, 4, 6, 8, 10, 12, 14, 16};

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

    std::ofstream outFile("threads_vs_time.csv");
    outFile << "ThreadCount,ExecutionTime\n";

    for (int n_threads : THREAD_COUNTS) {
        omp_set_num_threads(n_threads);

        int totalParams = NUM_QUBITS * 2 * NUM_LAYERS;
        std::vector<qreal> angles = generateRandomAngles(totalParams);

        double start = omp_get_wtime();

        #pragma omp parallel for
        for (int rep = 0; rep < 5; ++rep) {
            Qureg qureg = createDensityQureg(NUM_QUBITS);
            initZeroState(qureg);

            for (int layer = 0; layer < NUM_LAYERS; ++layer) {
                applyHEALayer(qureg, angles, layer, FIXED_NOISE);
            }

            destroyQureg(qureg);
        }

        double end = omp_get_wtime();
        double avg_time = (end - start) / 5.0;

        outFile << n_threads << "," << avg_time << "\n";
    }

    outFile.close();
    printf("Threads vs Time simulation complete.\n");

    return 0;
}
