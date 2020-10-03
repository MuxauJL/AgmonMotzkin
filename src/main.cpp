#include <iostream>
#include <numeric>

#include "AgmonMotzkin.h"
#include <Timer.h>

void appendObjectiveToConstraints(const ClpSimplex& model, AgmonMotzkin& rAgmonMotzkin,
    double dLowerBound, double dUpperBound)
{
    const double* pObjCoefficients = model.getObjCoefficients();
    int iNumberColumns = model.getNumCols();
    std::vector<int> indices(iNumberColumns);
    std::iota(indices.begin(), indices.end(), 0);
    rAgmonMotzkin.appendConstraint(iNumberColumns, indices.data(), pObjCoefficients,
        dLowerBound, dUpperBound);
}

int main()
{
    ClpSimplex model;
    int status = model.readMps("C:/Users/Mikhail/Desktop/repo/LinearProgramming/data/netlib/afiro.mps");
    //int status = model.readMps("C:/Users/Mikhail/Desktop/repo/LinearProgramming/data/netlib/OSA-60.mps");
    if (status)
    {
        return status;
    }

    {
        Timer t("Clp primal time: ");
        model.primal();
    }


    AgmonMotzkin agmonMotzkin;
    if (agmonMotzkin.init(model) == false)
    {
        return -1;
    }

    // optimal for afiro.mps
    //appendObjectiveToConstraints(model, agmonMotzkin, -464.77, -464.75);

    {
        Timer t("Agmon-Motzkin algorithm time: ");
        while (agmonMotzkin.nextApproximation(1e-8));
    }
    double dObj = agmonMotzkin.getCurrentObjectiveValue();
    printf("Agmon-Motzkin algorithm objective value = %f", dObj);

    return 0;
}
