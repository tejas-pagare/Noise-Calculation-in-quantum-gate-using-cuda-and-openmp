/** @file
 * API signatures for effecting decohering channels upon Quregs
 * which are instantiated as density matrices.
 * 
 * @author Tyson Jones
 * 
 * @defgroup decoherence Decoherence
 * @ingroup api
 * @brief Functions for effecting decoherence channels upon density matrices.
 * @{
 */

#ifndef DECOHERENCE_H
#define DECOHERENCE_H

#include "types.h"
#include "qureg.h"
#include "channels.h"


/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif

/// @notdoced
/// @notvalidated
void mixDephasing(Qureg qureg, int qubit, qreal prob);

/// @notdoced
/// @notvalidated
void mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal prob);

/// @notdoced
/// @notvalidated
void mixDepolarising(Qureg qureg, int qubit, qreal prob);

/// @notdoced
/// @notvalidated
void mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal prob);

/// @notdoced
/// @notvalidated
void mixDamping(Qureg qureg, int qubit, qreal prob);

/// @notdoced
/// @notvalidated
void mixPaulis(Qureg qureg, int qubit, qreal probX, qreal probY, qreal probZ);

/// @notdoced
/// @notvalidated
void mixQureg(Qureg qureg, Qureg other, qreal prob);

/// @notdoced
/// @notvalidated
void mixKrausMap(Qureg qureg, int* targets, int numTargets, KrausMap map);

/// @notdoced
/// @notvalidated
void mixSuperOp(Qureg qureg, int* targets, int numTargets, SuperOp superop);


void mixDephasingTracked(Qureg qureg, int qubit, qreal prob, qreal* errorOut) ;
// void mixTwoQubitDephasingTracked(Qureg qureg, int qubit1, int qubit2, qreal prob, qreal* errorOut);
// void mixDepolarisingTracked(Qureg qureg, int qubit, qreal prob, qreal* errorOut);
// void mixTwoQubitDepolarisingTracked(Qureg qureg, int qubit1, int qubit2, qreal prob, qreal* errorOut) ;
// void mixDampingTracked(Qureg qureg, int qubit, qreal prob, qreal* errorOut) ;

// end de-mangler
#ifdef __cplusplus
}
#endif



/*
 * C++ OVERLOADS
 *
 * which are only accessible to C++ binaries, and accept
 * arguments more natural to C++ (e.g. std::vector). These
 * are included in the file-wide doxygen group (no subgroups).
 */

#ifdef __cplusplus

#include <vector>

/// @nottested
/// @notdoced
/// @notvalidated
/// @cpponly
void mixKrausMap(Qureg qureg, std::vector<int> targets, KrausMap map);

/// @nottested
/// @notdoced
/// @notvalidated
/// @cpponly
void mixSuperOp(Qureg qureg, std::vector<int> targets, SuperOp superop);

#endif // __cplusplus



#endif // DECOHERENCE_H

/** @} */ // (end file-wide doxygen defgroup)
