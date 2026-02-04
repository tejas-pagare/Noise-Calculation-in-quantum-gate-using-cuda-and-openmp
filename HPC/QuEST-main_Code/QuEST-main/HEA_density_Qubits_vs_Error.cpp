#include "quest.h"
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>

using qreal = qreal;

const int NUM_LAYERS_PER_QUBIT = 5; // Scale layers linearly
const std::vector<int> QUBIT_COUNTS = {4, 6, 8, 10, 12, 14, 16, 18, 20};

std::vector<qreal> generateRandomAngles(int numParams) {
    std::vector<qreal> angles(numParams);
    for (int i = 0; i < numParams; ++i)
        angles[i] = ((qreal)rand() / RAND_MAX) * M_PI;
    return angles;
}

void applyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex, qreal noise) {
    int offset = layerIndex * qureg.numQubits * 2;

    for (int qubit = 0; qubit < qureg.numQubits; ++qubit) {
        applyRotateX(qureg, qubit, angles[offset + qubit * 2]);
        applyRotateY(qureg, qubit, angles[offset + qubit * 2 + 1]);
        qreal error;
        mixDephasing(qureg, qubit, noise, &error);
    }

    for (int qubit = 0; qubit < qureg.numQubits - 1; ++qubit) {
        applyControlledMultiQubitNot(qureg, qubit, new int[1]{qubit + 1}, 1);
    }
}

int main() {
    srand(time(0));
    initQuESTEnv();

    std::ofstream outFile("qubits_vs_error.csv");
    outFile << "NumQubits,NumLayers,AverageError\n";

    #pragma omp parallel for schedule(dynamic)
    for (size_t idx = 0; idx < QUBIT_COUNTS.size(); ++idx) {
        int numQubits = QUBIT_COUNTS[idx];
        int numLayers = numQubits * NUM_LAYERS_PER_QUBIT;

        Qureg qureg = createDensityQureg(numQubits);
        int totalParams = numQubits * 2 * numLayers;
        std::vector<qreal> angles = generateRandomAngles(totalParams);

        initZeroState(qureg);
        for (int layer = 0; layer < numLayers; ++layer) {
            applyHEALayer(qureg, angles, layer, 0.10); // fixed noise
        }

        Qureg target = createDensityQureg(numQubits);
        initZeroState(target);
        qreal fidelity = calcFidelity(qureg, target);
        qreal error = 1.0 - fidelity;

        outFile << numQubits << "," << numLayers << "," << error << "\n";
        destroyQureg(qureg);
        destroyQureg(target);
    }

    outFile.close();
    printf("Qubits vs Error simulation complete.\n");

    return 0;
}
