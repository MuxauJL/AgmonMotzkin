#ifndef MV_AgmonMotzkin_H_

#include <vector>
#include <unordered_map>
#include <ClpSimplex.hpp>
#include <CoinPackedVector.hpp>

class AgmonMotzkin
{
public:
    AgmonMotzkin() = default;
    ~AgmonMotzkin() = default;
    bool init(const ClpSimplex& rModel);
    bool setCurrentPoint(const std::vector<double>& rCurrentPoint);
    const std::vector<double>& getCurrentPoint();
    double getCurrentObjectiveValue();
    void appendConstraint(int iNumCoefficients, const int* pIndices, const double* pCoefficients, double dLowerBound, double dUpperBound = COIN_DBL_MAX);
    void printMatrix();
    std::pair<int, double> prepareStep();
    std::pair<int, double> prepareStep(int iProcessId, int iNumProcesses);
    void updateCurrentPoint(int iRow, double dStep);
private:
    void convertToStandardFormMin(std::vector<double>& rRowUppers);
    double coefficientSquaresSum(int iRowIndex);

private:
    CoinPackedMatrix matrix_;
    int optimizationDirection_ = 0; // 0 - not set, 1 - min, -1 - max
    std::vector<double> objCoefficients_;
    // rowUppers_ and colLowers_ are equal to +inf
    // since we store the matrix_ of the minimization problem in a standard form.
    std::vector<double> rowLowers_;
    std::vector<double> colLowers_;
    std::vector<double> x_;
    std::unordered_map<int, double> coefficientSquaresSums_; // cache
};

#endif // MV_AgmonMotzkin_H_
