#include <iostream>
#include "AgmonMotzkin.h"
#include <numeric>

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
    if (status)
    {
        return status;
    }
    model.primal();


    AgmonMotzkin agmonMotzkin;
    if (agmonMotzkin.init(model) == false)
    {
        return -1;
    }

    appendObjectiveToConstraints(model, agmonMotzkin, -464.77, -464.75);

    while (agmonMotzkin.nextApproximation(1e-8));
    auto x = agmonMotzkin.getCurrentPoint();
    double dObj = agmonMotzkin.getCurrentObjectiveValue();

    const double* columnPrimal = model.getColSolution();
    int iNumberColumns = model.getNumCols();
    for (int iCol = 0; iCol < iNumberColumns; ++iCol)
    {
        printf("%f %f\n", x[iCol], columnPrimal[iCol]);
    }
    printf("Agmon-Motzkin algorithm objective value = %f", dObj);

    return 0;
}
