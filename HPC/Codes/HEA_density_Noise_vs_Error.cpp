#include "quest.h"
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>

using qreal = qreal;

const int NUM_QUBITS = 10;
const int NUM_LAYERS = 10;

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

    // Single-qubit rotations
    for (int qubit = 0; qubit < NUM_QUBITS; ++qubit) {
        applyRotateX(qureg, qubit, angles[offset + qubit * 2]);
        applyRotateY(qureg, qubit, angles[offset + qubit * 2 + 1]);

        qreal error;
        mixDephasing(qureg, qubit, noise); // Apply noise
    }

    // Entanglement layer (linear chain of CNOTs)
    for (int qubit = 0; qubit < NUM_QUBITS - 1; ++qubit) {
        applyControlledMultiQubitNot(qureg, qubit, new int[1]{qubit + 1}, 1);
    }
}

int main() {
    srand(time(0));
    initQuESTEnv();

    std::ofstream outFile("noise_vs_error_density.csv");
    outFile << "NoiseLevel,AverageError\n";

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < NOISE_LEVELS.size(); ++i) {
        qreal noise = NOISE_LEVELS[i];
        Qureg qureg = createDensityQureg(NUM_QUBITS);

        initZeroState(qureg);
        int totalParams = NUM_QUBITS * 2 * NUM_LAYERS;
        std::vector<qreal> angles = generateRandomAngles(totalParams);

        for (int layer = 0; layer < NUM_LAYERS; ++layer) {
            applyHEALayer(qureg, angles, layer, noise);
        }

        // Create a pure-state register for the ideal state
Qureg target_pure = createQureg(NUM_QUBITS);
initZeroState(target_pure);

// Fidelity between density matrix and pure state IS supported
qreal fidelity = calcFidelity(qureg, target_pure);
qreal error = 1.0 - fidelity;


        outFile << noise << "," << error << "\n";

        destroyQureg(qureg);
        destroyQureg(target_pure);
    }

    outFile.close();
    printf("Noise vs Error simulation complete.\n");

    return 0;
}
