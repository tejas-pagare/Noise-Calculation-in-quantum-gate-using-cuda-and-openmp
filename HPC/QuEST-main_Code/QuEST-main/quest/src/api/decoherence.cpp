/** @file
 * API definitions for effecting decohering channels upon Quregs
 * which are instantiated as density matrices.
 * 
 * @author Tyson Jones
 * @author Balint Koczor (prototyped v3 mixKrausMap)
 * @author Nicolas Vogt (prototyped v3 mixDamping)
 */

#include "../../include/types.h"




#include "../../include/qureg.h"




#include "../../include/channels.h"





#include "../src/core/validation.hpp"

#include "../src/core/localiser.hpp"

#include "../src/core/utilities.hpp"


#include <vector>
using std::vector;



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
extern "C" {


void mixDephasing(Qureg qureg, int qubit, qreal prob) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_oneQubitDepashingProb(prob, __func__);

    // permit but do not change non-decohering statevecs
    if (prob == 0) 
        return;
    
    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_oneQubitDephasing(qureg, qubit, prob);
}


void mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal prob) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, qubit1, qubit2, __func__);
    validate_twoQubitDepashingProb(prob, __func__);

    // permit but do not change non-decohering statevecs
    if (prob == 0) 
        return;
    
    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_twoQubitDephasing(qureg, qubit1, qubit2, prob);
}


void mixDepolarising(Qureg qureg, int qubit, qreal prob) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_oneQubitDepolarisingProb(prob, __func__);

    // permit but do not change non-decohering statevecs
    if (prob == 0) 
        return;
    
    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_oneQubitDepolarising(qureg, qubit, prob);
}


void mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal prob) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, qubit1, qubit2, __func__);
    validate_twoQubitDepolarisingProb(prob, __func__);

    // permit but do not change non-decohering statevecs
    if (prob == 0) 
        return;
    
    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_twoQubitDepolarising(qureg, qubit1, qubit2, prob);
}


void mixDamping(Qureg qureg, int qubit, qreal prob) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_oneQubitDampingProb(prob, __func__);

    // permit but do not change non-decohering statevecs
    if (prob == 0) 
        return;
    
    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_oneQubitDamping(qureg, qubit, prob);
}


void mixPaulis(Qureg qureg, int qubit, qreal probX, qreal probY, qreal probZ) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_oneQubitPauliChannelProbs(probX, probY, probZ, __func__);

    // permit but do not change non-decohering statevecs
    if (probX == 0 && probY == 0 && probZ == 0)
        return;

    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_oneQubitPauliChannel(qureg, qubit, probX, probY, probZ);
}


void mixKrausMap(Qureg qureg, int* qubits, int numQubits, KrausMap map) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, qubits, numQubits, __func__);
    validate_mixedAmpsFitInNode(qureg, 2*numQubits, __func__); // superop acts on 2x
    validate_krausMapIsCPTP(map, __func__); // also checks fields and is-sync
    validate_krausMapMatchesTargets(map, numQubits, __func__);

    localiser_densmatr_krausMap(qureg, map, util_getVector(qubits, numQubits));
}


void mixQureg(Qureg outQureg, Qureg inQureg, qreal inProb) {
    validate_quregFields(outQureg, __func__);
    validate_quregFields(inQureg, __func__);
    validate_probability(inProb, __func__);
    validate_quregsCanBeMixed(outQureg, inQureg, __func__); // checks outQureg is densmatr

    qreal outProb = 1 - inProb;
    localiser_densmatr_mixQureg(outProb, outQureg, inProb, inQureg);
}


void mixSuperOp(Qureg qureg, int* targets, int numTargets, SuperOp superop) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_superOpFields(superop, __func__);
    validate_superOpIsSynced(superop, __func__);
    validate_superOpDimMatchesTargs(superop, numTargets, __func__);
    validate_mixedAmpsFitInNode(qureg, 2*numTargets, __func__); // superop acts on 2x

    localiser_densmatr_superoperator(qureg, superop, util_getVector(targets, numTargets));
}

// void mixDephasingTracked(Qureg qureg, int qubit, qreal prob, qreal* errorOut) {
//     validate_quregFields(qureg, __func__);
//     validate_target(qureg, qubit, __func__);
//     validate_oneQubitDepashingProb(prob, __func__);

//     if (prob == 0) {
//         if (errorOut) *errorOut = 0;
//         return;
//     }

//     validate_quregIsDensityMatrix(qureg, __func__);

//     // Track error with parallelism
//     qreal localError = 0;
//     #pragma omp parallel for reduction(+:localError)
//     for (long long int i = 0; i < qureg.numAmps; i++) {
//         qreal ampBefore = abs(qureg.rhoRe[i][i]);
//         // localiser applies the decoherence
//         localiser_densmatr_oneQubitDephasing(qureg, qubit, prob);
//         qreal ampAfter = abs(qureg.rhoRe[i][i]);
//         localError += fabs(ampBefore - ampAfter);
//     }

//     if (errorOut)
//         *errorOut = localError;
// }
// void mixDephasingTracked(Qureg qureg, int qubit, qreal prob, qreal* errorOut) {
//     validate_quregFields(qureg, __func__);
//     validate_target(qureg, qubit, __func__);
//     validate_oneQubitDepashingProb(prob, __func__);

//     if (prob == 0) {
//         if (errorOut) *errorOut = 0;
//         return;
//     }

//     validate_quregIsDensityMatrix(qureg, __func__);

//     long long int dim = qureg.numAmpsPerNode;  // total number of elements in rho
//     long long int N = (long long int) sqrt(dim);

//     // Capture diagonal values BEFORE decoherence
//     qreal* diagBefore = (qreal*) malloc(N * sizeof(qreal));
//     #pragma omp parallel for
//     for (long long int i = 0; i < N; i++) {
//         long long int diagIdx = i*N + i;
//         diagBefore[i] = fabs(qureg.cpuAmps[diagIdx].real);  // only real part
//     }

//     // Apply decoherence
//     localiser_densmatr_oneQubitDephasing(qureg, qubit, prob);

//     // Track error: compare AFTER vs BEFORE
//     qreal localError = 0;
//     #pragma omp parallel for reduction(+:localError)
//     for (long long int i = 0; i < N; i++) {
//         long long int diagIdx = i*N + i;
//         qreal diagAfter = fabs(qureg.cpuAmps[diagIdx].real);
//         localError += fabs(diagAfter - diagBefore[i]);
//     }

//     if (errorOut)
//         *errorOut = localError;

//     free(diagBefore);
// }
// void mixDephasingTracked(Qureg qureg, int qubit, qreal prob, qreal* errorOut) {
    
//     validate_quregFields(qureg, __func__);
//     validate_target(qureg, qubit, __func__);
//     validate_oneQubitDephasingProb(prob, __func__);
    
//     if (prob == 0) {
//         if (errorOut) *errorOut = 0;
//         return;
//     }

//     validate_quregIsDensityMatrix(qureg, __func__);

//     long long int numQubits = qureg.numQubitsRepresented;
//     long long int numAmps = qureg.numAmpsTotal;
//     long long int numCols = qureg.numAmpsPerRow;

//     // Allocate array to store real part of diagonal elements before dephasing
//     qreal* diagBefore = (qreal*) malloc(numAmps * sizeof(qreal));
//     if (diagBefore == NULL) {
//         fprintf(stderr, "ERROR: mixDephasingTracked failed to allocate memory\n");
//         exit(1);
//     }

//     // Store only diagonal (i == j) elements
//     for (long long int i = 0; i < numAmps; i++) {
//         long long int diagIdx = i * numCols + i;  // row * stride + col
//         diagBefore[i] = qureg.cpuAmps[diagIdx].real();
//     }

//     // Apply standard dephasing
//     mixDephasing(qureg, qubit, prob);

//     // Compute total absolute error in diagonals
//     // qreal totalError = 0;
//     qreal localError = 0;
//     #pragma omp parallel for reduction(+:localError)
//     for (long long int i = 0; i < numAmps; i++) {
//         long long int diagIdx = i * numCols + i;
//         qreal diagAfter = qureg.cpuAmps[diagIdx].real();
//         localError += fabs(diagBefore[i] - diagAfter);
//     }

//     free(diagBefore);

//     if (errorOut)
//         *errorOut = localError;
// }
// void mixDephasingTracked(Qureg qureg, int qubit, qreal prob, qreal* errorOut) {
//     validate_quregFields(qureg, __func__);
//     validate_target(qureg, qubit, __func__);
//     validate_oneQubitDepashingProb(prob, __func__);

//     if (prob == 0) {
//         if (errorOut) *errorOut = 0;
//         return;
//     }

//     validate_quregIsDensityMatrix(qureg, __func__);

//     long long int dim = qureg.numAmpsPerNode;  // total number of elements in rho
//     long long int N = (long long int) sqrt(dim);

//     // Capture diagonal values BEFORE decoherence
//     qreal* diagBefore = (qreal*) malloc(N * sizeof(qreal));
//     #pragma omp parallel for
//     for (long long int i = 0; i < N; i++) {
//         long long int diagIdx = i*N + i;
//         diagBefore[i] = fabs(qureg.cpuAmps[diagIdx].real());  // only real part
//     }

//     // Apply decoherence
//     localiser_densmatr_oneQubitDephasing(qureg, qubit, prob);

//     // Track error: compare AFTER vs BEFORE
//     qreal localError = 0;
//     #pragma omp parallel for reduction(+:localError)
//     for (long long int i = 0; i < N; i++) {
//         long long int diagIdx = i*N + i;
//         qreal diagAfter = fabs(qureg.cpuAmps[diagIdx].real());
//         localError += fabs(diagAfter - diagBefore[i]);
//     }

//     if (errorOut)
//         *errorOut = localError;

//     free(diagBefore);
// }

void mixDephasingTracked(Qureg qureg, int qubit, qreal prob, qreal* errorOut) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_oneQubitDepashingProb(prob, __func__);

    if (prob == 0) {
        if (errorOut) *errorOut = 0;
        return;
    }

    validate_quregIsDensityMatrix(qureg, __func__);

    long long int dim = qureg.numAmpsPerNode;  // total number of elements in rho

    // Capture full real part of rho BEFORE
    qreal* beforeRe = (qreal*) malloc(dim * sizeof(qreal));
    #pragma omp parallel for
    for (long long int i = 0; i < dim; i++) {
        beforeRe[i] = qureg.cpuAmps[i].real();
    }

    // Apply decoherence
    localiser_densmatr_oneQubitDephasing(qureg, qubit, prob);

    // Compare all off-diagonal real parts
    long long int N = (long long int) sqrt(dim);
    qreal localError = 0;
    #pragma omp parallel for reduction(+:localError)
    for (long long int i = 0; i < N; i++) {
        for (long long int j = 0; j < N; j++) {
            if (i == j) continue; // skip diagonal

            long long int idx = i*N + j;
            qreal afterRe = qureg.cpuAmps[idx].real();
            localError += fabs(afterRe - beforeRe[idx]);
        }
    }

    if (errorOut)
        *errorOut = localError;

    free(beforeRe);
}

// void mixTwoQubitDephasingTracked(Qureg qureg, int qubit1, int qubit2, qreal prob, qreal* errorOut) {
//     validate_quregFields(qureg, __func__);
//     validate_twoTargets(qureg, qubit1, qubit2, __func__);
//     validate_twoQubitDepashingProb(prob, __func__);

//     // permit but do not change non-decohering statevecs
//      if (prob == 0) {
//         if (errorOut) *errorOut = 0;
//         return;
//     }

//     validate_quregIsDensityMatrix(qureg, __func__);

//     // Track error with parallelism
//     qreal localError = 0;
//     #pragma omp parallel for reduction(+:localError)
//     for (long long int i = 0; i < qureg.numAmps; i++) {
//         qreal ampBefore = abs(qureg.rhoRe[i][i]);
//         // localiser applies the decoherence
//         localiser_densmatr_twoQubitDephasing(qureg, qubit1, qubit2, prob);qreal ampAfter = abs(qureg.rhoRe[i][i]);
//         localError += fabs(ampBefore - ampAfter);
//     }

//     if (errorOut)
//         *errorOut = localError;
// }


// void mixDepolarisingTracked(Qureg qureg, int qubit, qreal prob, qreal* errorOut) {
//     validate_quregFields(qureg, __func__);
//     validate_target(qureg, qubit, __func__);
//     validate_oneQubitDepolarisingProb(prob, __func__);

//     // permit but do not change non-decohering statevecs
//     if (prob == 0) {
//         if (errorOut) *errorOut = 0;
//         return;
//     }

//     validate_quregIsDensityMatrix(qureg, __func__);

//     // Track error with parallelism
//     qreal localError = 0;
//     #pragma omp parallel for reduction(+:localError)
//     for (long long int i = 0; i < qureg.numAmps; i++) {
//         qreal ampBefore = abs(qureg.rhoRe[i][i]);
//         // localiser applies the decoherence
//         localiser_densmatr_oneQubitDepolarising(qureg, qubit, prob); qreal ampAfter = abs(qureg.rhoRe[i][i]);
//         localError += fabs(ampBefore - ampAfter);
//     }

//     if (errorOut)
//         *errorOut = localError;
// }

// void mixTwoQubitDepolarisingTracked(Qureg qureg, int qubit1, int qubit2, qreal prob, qreal* errorOut) {
//     validate_quregFields(qureg, __func__);
//     validate_twoTargets(qureg, qubit1, qubit2, __func__);
//     validate_twoQubitDepolarisingProb(prob, __func__);

//     // permit but do not change non-decohering statevecs
//     if (prob == 0) {
//         if (errorOut) *errorOut = 0;
//         return;
//     }

//     validate_quregIsDensityMatrix(qureg, __func__);

//     // Track error with parallelism
//     qreal localError = 0;
//     #pragma omp parallel for reduction(+:localError)
//     for (long long int i = 0; i < qureg.numAmps; i++) {
//         qreal ampBefore = abs(qureg.rhoRe[i][i]);
//         // localiser applies the decoherence
//         localiser_densmatr_twoQubitDepolarising(qureg, qubit1, qubit2, prob);
// qreal ampAfter = abs(qureg.rhoRe[i][i]);
//         localError += fabs(ampBefore - ampAfter);
//     }

//     if (errorOut)
//         *errorOut = localError;
// }

// void mixDampingTracked(Qureg qureg, int qubit, qreal prob, qreal* errorOut) {
//     validate_quregFields(qureg, __func__);
//     validate_target(qureg, qubit, __func__);
//     validate_oneQubitDampingProb(prob, __func__);

//     // permit but do not change non-decohering statevecs
//     if (prob == 0) {
//         if (errorOut) *errorOut = 0;
//         return;
//     }

//     validate_quregIsDensityMatrix(qureg, __func__);

//     // Track error with parallelism
//     qreal localError = 0;
//     #pragma omp parallel for reduction(+:localError)
//     for (long long int i = 0; i < qureg.numAmps; i++) {
//         qreal ampBefore = abs(qureg.rhoRe[i][i]);
//         // localiser applies the decoherence
//         localiser_densmatr_oneQubitDamping(qureg, qubit, prob);
// qreal ampAfter = abs(qureg.rhoRe[i][i]);
//         localError += fabs(ampBefore - ampAfter);
//     }

//     if (errorOut)
//         *errorOut = localError;
// }

} // end de-mangler



/*
 * C++ OVERLOADS
 *
 * which are only accessible to C++ binaries, and accept
 * arguments more natural to C++ (e.g. std::vector). We
 * manually add these to their respective Doxygen doc groups.
 */

#ifdef __cplusplus

void mixKrausMap(Qureg qureg, vector<int> targets, KrausMap map) {
    mixKrausMap(qureg, targets.data(), targets.size(), map);
}

void mixSuperOp(Qureg qureg, vector<int> targets, SuperOp superop) {
    mixSuperOp(qureg, targets.data(), targets.size(), superop);
}

#endif // __cplusplus
