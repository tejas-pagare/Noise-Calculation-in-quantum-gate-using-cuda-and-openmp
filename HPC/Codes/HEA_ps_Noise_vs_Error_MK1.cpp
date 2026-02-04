// noise_vs_error.cpp
#include "quest.h"
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>

using qreal = qreal;

const int NUM_QUBITS = 20;
const int NUM_LAYERS = 50;

std::vector<qreal> NOISE_LEVELS = {
    0.00, 0.02, 0.04, 0.06,
    0.08, 0.10, 0.12, 0.14,
    0.16, 0.18, 0.20, 0.22,
    0.24, 0.26, 0.28, 0.30,
    0.32, 0.34, 0.36, 0.38,
    0.40, 0.42, 0.44, 0.46,
    0.48, 0.50
};

std::vector<qreal> generateRandomAngles(int numParams) {
    std::vector<qreal> angles(numParams);
    for (int i = 0; i < numParams; ++i)
        angles[i] = ((qreal)rand() / RAND_MAX) * M_PI;
    return angles;
}

void applyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex, qreal noise) {
    int offset = layerIndex * NUM_QUBITS * 2;

    for (int qubit = 0; qubit < NUM_QUBITS; ++qubit) {
        qreal rx_angle = angles[offset + qubit * 2];
        qreal ry_angle = angles[offset + qubit * 2 + 1];

        rx_angle *= (1 + noise * ((qreal)rand() / RAND_MAX - 0.5));
        ry_angle *= (1 + noise * ((qreal)rand() / RAND_MAX - 0.5));

        applyRotateX(qureg, qubit, rx_angle);
        applyRotateY(qureg, qubit, ry_angle);
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

    std::ofstream outFile("noise_vs_error.csv");
    outFile << "NoiseLevel,AverageError\n";

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < NOISE_LEVELS.size(); ++i) {
        qreal noise = NOISE_LEVELS[i];
        qreal totalError = 0.0;
        const int NUM_RUNS = 5;

        for (int run = 0; run < NUM_RUNS; ++run) {
            Qureg qureg = createQureg(NUM_QUBITS);
            initZeroState(qureg);

            for (int layer = 0; layer < NUM_LAYERS; ++layer) {
                applyHEALayer(qureg, fixedAngles, layer, noise);
            }

            Qureg target = createQureg(NUM_QUBITS);
            initZeroState(target);

            qreal fidelity = calcFidelity(qureg, target);
            qreal error = 1.0 - fidelity;
            totalError += error;

            destroyQureg(qureg);
            destroyQureg(target);
        }

        qreal avgError = totalError / NUM_RUNS;

        // Ensure thread-safe write to output file
        #pragma omp critical
        {
            outFile << noise << "," << avgError << "\n";
        }
    }

    outFile.close();
    printf("Noise vs Error simulation complete.\n");
    return 0;
}
