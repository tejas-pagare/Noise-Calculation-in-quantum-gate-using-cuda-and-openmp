#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include "quest.h"

// ----- Parameters -----
const int n_qubits = 4;
const int min_depth = 1;
const int max_depth = 5;
const double EPS = 1e-3;       // Finite difference step for gradient
const std::vector<double> noise_levels = {0.0, 0.01, 0.05, 0.1};

std::default_random_engine rng(std::random_device{}());
std::normal_distribution<double> gauss(0.0, 1.0);

// ----- Helper: Add Gaussian noise to parameter -----
double addGaussianNoise(double param, double noise_level) {
    return param + noise_level * gauss(rng);
}

// ----- Build HEA Circuit: rotations + entangling CX (linear chain) -----
void buildHEACircuit(Qureg qureg, const std::vector<double>& angles, int depth) {
    int idx = 0;
    for (int layer = 0; layer < depth; layer++) {
        // Single-qubit rotations (Ry)
        for (int q = 0; q < n_qubits; q++) {
            rotateY(qureg, q, angles[idx++]);
        }
        // Entangling CX gates in a chain
        for (int q = 0; q < n_qubits - 1; q++) {
            controlledNot(qureg, q, q + 1);
        }
    }
}

// ----- Cost function: expectation value of Z on qubit 0 -----
double costFunction(const std::vector<double>& angles, int depth, double noise_level) {
    // Initialize state |0...0>
    Qureg qureg = createQureg(n_qubits, globalQubitRegister());
    initZeroState(qureg);

    // Apply noisy angles
    std::vector<double> noisy_angles(angles.size());
    for (size_t i = 0; i < angles.size(); i++) {
        noisy_angles[i] = addGaussianNoise(angles[i], noise_level);
    }

    buildHEACircuit(qureg, noisy_angles, depth);

    // Measure <Z_0> = Probability(0) - Probability(1) on qubit 0
    double p0 = getProbAmp(qureg, 0);     // Amplitude squared of |0...0> (binary state 0)
    double p1 = 1.0 - p0;                  // Probability qubit 0 = 1 is 1 - p0 (approximately)
    double expZ0 = p0 - p1;

    destroyQureg(qureg, globalQubitRegister());
    return expZ0;
}

// ----- Compute all gradients via finite differences -----
std::vector<double> computeGradients(const std::vector<double>& angles, int depth, double noise_level) {
    int n_params = angles.size();
    std::vector<double> grads(n_params);

    for (int i = 0; i < n_params; i++) {
        std::vector<double> angles_plus = angles;
        std::vector<double> angles_minus = angles;

        angles_plus[i] += EPS;
        angles_minus[i] -= EPS;

        double c_plus = costFunction(angles_plus, depth, noise_level);
        double c_minus = costFunction(angles_minus, depth, noise_level);
        grads[i] = (c_plus - c_minus) / (2 * EPS);
    }
    return grads;
}

// ----- Generate initial angles (small random values near zero) -----
std::vector<double> generateInitialAngles(int n_qubits, int depth) {
    int n_params = n_qubits * depth;
    std::vector<double> angles(n_params);
    std::normal_distribution<double> smallNoise(0.0, 0.01);
    for (int i = 0; i < n_params; i++) {
        angles[i] = smallNoise(rng);
    }
    return angles;
}

// ----- Main benchmarking loop -----
int main(int argc, char *argv[]) {
    // Initialize QuEST environment
    QuESTEnv env = createQuESTEnv();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Depth,NoiseLevel,Cost,GradNorm\n";

    for (double noise_level : noise_levels) {
        for (int depth = min_depth; depth <= max_depth; depth++) {
            // Generate initial parameters
            std::vector<double> angles = generateInitialAngles(n_qubits, depth);

            // Evaluate cost function
            double cost = costFunction(angles, depth, noise_level);

            // Compute gradients
            std::vector<double> gradients = computeGradients(angles, depth, noise_level);

            // Compute gradient norm
            double grad_norm = 0.0;
            for (double g : gradients) grad_norm += g * g;
            grad_norm = std::sqrt(grad_norm);

            // Output results as CSV format
            std::cout << depth << "," << noise_level << "," << cost << "," << grad_norm << "\n";
        }
    }

    destroyQuESTEnv(env);

    return 0;
}

