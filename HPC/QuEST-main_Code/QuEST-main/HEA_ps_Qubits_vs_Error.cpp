// qubits_vs_error_scaled.cpp
#include "quest.h"
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>

using qreal = qreal;

const qreal FIXED_NOISE = 0.05;

std::vector<std::pair<int, int>> QUBIT_LAYER_COMBOS = {
    {4, 10}, {6, 15}, {8, 20}, {10, 25},
    {12, 30}, {14, 35}, {16, 40}, {18, 45}, {20, 50}, {22, 55}
};

std::vector<qreal> generateRandomAngles(int numParams) {
    std::vector<qreal> angles(numParams);
    for (int i = 0; i < numParams; ++i)
        angles[i] = ((qreal)rand() / RAND_MAX) * M_PI;
    return angles;
}

void applyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex, qreal noise) {
    int offset = layerIndex * qureg.numQubits * 2;

    for (int qubit = 0; qubit < qureg.numQubits; ++qubit) {
        qreal rx_angle = angles[offset + qubit * 2];
        qreal ry_angle = angles[offset + qubit * 2 + 1];

        rx_angle *= (1 + noise * ((qreal)rand() / RAND_MAX - 0.5));
        ry_angle *= (1 + noise * ((qreal)rand() / RAND_MAX - 0.5));

        applyRotateX(qureg, qubit, rx_angle);
        applyRotateY(qureg, qubit, ry_angle);
    }

    for (int qubit = 0; qubit < qureg.numQubits - 1; ++qubit) {
        applyControlledMultiQubitNot(qureg, qubit, new int[1]{qubit + 1}, 1);
    }
}

int main() {
    srand(time(0));
    initQuESTEnv();

    std::ofstream outFile("qubits_vs_error_scaled.csv");
    outFile << "NumQubits,NumLayers,AverageError\n";

    #pragma omp parallel for schedule(dynamic)
    for (size_t idx = 0; idx < QUBIT_LAYER_COMBOS.size(); ++idx) {
        int numQubits = QUBIT_LAYER_COMBOS[idx].first;
        int numLayers = QUBIT_LAYER_COMBOS[idx].second;

        int totalParams = numQubits * 2 * numLayers;
        std::vector<qreal> fixedAngles = generateRandomAngles(totalParams);

        Qureg qureg = createQureg(numQubits);

        initZeroState(qureg);
        for (int layer = 0; layer < numLayers; ++layer) {
            applyHEALayer(qureg, fixedAngles, layer, FIXED_NOISE);
        }

       Qureg target = createQureg(numQubits); // or 'numQubits' if variable
initZeroState(target);
qreal fidelity = calcFidelity(qureg, target);
qreal error = 1.0 - fidelity;
destroyQureg(target);

        outFile << numQubits << "," << numLayers << "," << error << "\n";
        destroyQureg(qureg);
    }

    printf("Scaled Qubits vs Error simulation complete.\n");
    return 0;
}
