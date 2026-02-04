#include "quest.h"
#include <fstream>

void simulateBell(Qureg qureg, qreal noise, bool useTracked, std::ofstream& outFile) {
    applyHadamard(qureg, 0);
    int target = 1;
applyControlledMultiQubitNot(qureg, 0, &target, 1);


    qreal error = 0;
    if (useTracked)
        mixDephasingTracked(qureg, 0, noise, &error);
    else
        mixDephasing(qureg, 0, noise);

    outFile << "Bell," << noise << "," << useTracked << "," << error << "\n";
}

int main() {
    initQuESTEnv() ;
    QuESTEnv env = getQuESTEnv();
    Qureg qureg = createDensityQureg(2);
    initZeroState(qureg);

    std::ofstream out("results.csv", std::ios::app);
    simulateBell(qureg, 0.3, false, out); // untracked
    simulateBell(qureg, 0.3, true, out);  // tracked
    destroyQureg(qureg);
    // destroyQuESTEnv(env);
}
