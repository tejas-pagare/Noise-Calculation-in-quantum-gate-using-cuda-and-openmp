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
    0.00, 0.10, 0.20, 0.30, 0.40, 0.50
};

std::vector<qreal> generateRandomAngles(int numParams) {
    std::vector<qreal> angles(numParams);
    for (int i = 0; i < numParams; ++i)
        angles[i] = ((qreal)rand() / RAND_MAX) * M_PI;
    return angles;
}

void applyHEALayer(Qureg qureg, const std::vector<qreal>& angles, int layerIndex, qreal noise) {
    int offset = layerIndex * NUM_QUBITS * 2;

    // Single-qubit rotations + noise
    for (int qubit = 0; qubit < NUM_QUBITS; ++qubit) {
        applyRotateX(qureg, qubit, angles[offset + qubit * 2]);
        applyRotateY(qureg, qubit, angles[offset + qubit * 2 + 1]);
        mixDephasing(qureg, qubit, noise); // Apply dephasing noise
    }

    // Entanglement layer (CNOTs)
    for (int qubit = 0; qubit < NUM_QUBITS - 1; ++qubit) {
        applyCNOT(qureg, qubit, qubit + 1);
    }
}

int main() {
    srand(time(0));

    QuESTEnv env = createQuESTEnv();

    std::ofstream outFile("fidelity_per_layer.csv");
    outFile << "NoiseLevel,Layer,Fidelity\n";

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < NOISE_LEVELS.size(); ++i) {
        qreal noise = NOISE_LEVELS[i];

        Qureg qureg = createDensityQureg(NUM_QUBITS, env);
        initZeroState(qureg);

        int totalParams = NUM_QUBITS * 2 * NUM_LAYERS;
        std::vector<qreal> angles = generateRandomAngles(totalParams);

        Qureg target_pure = createQureg(NUM_QUBITS, env);
        initZeroState(target_pure);

        for (int layer = 0; layer < NUM_LAYERS; ++layer) {
            applyHEALayer(qureg, angles, layer, noise);

            qreal fidelity = calcFidelity(qureg, target_pure);

            #pragma omp critical
            {
                outFile << noise << "," << layer + 1 << "," << fidelity << "\n";
            }
        }

        destroyQureg(qureg);
        destroyQureg(target_pure);
    }

    destroyQuESTEnv(env);
    outFile.close();
    printf("Fidelity per layer simulation complete.\n");

    return 0;
}
