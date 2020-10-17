#include <Windows.h>

#include <iostream>
#include <numeric>
#include <fstream>
#include <filesystem>

#include "AgmonMotzkin.h"
#include <Timer.h>

namespace fs = std::filesystem;

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

bool saveBufferToBinFile(const std::string& sPathToSave, const char* pBuffer, size_t nBufferSize)
{
    std::ofstream fs(sPathToSave, std::ios::out | std::ios::binary | std::ios::trunc);
    if (fs.is_open())
    {
        fs.write(pBuffer, nBufferSize);
        fs.close();
        return true;
    }
    else
    {
        return false;
    }
}

template<typename T>
bool saveToBinFile(const std::string& sPathToSave, const std::vector<T>& rInput)
{
    return saveBufferToBinFile(sPathToSave,
        reinterpret_cast<const char*>(rInput.data()), rInput.size() * sizeof(T));
}

template<typename T>
bool loadBinFile(const std::string& sPathToFile, std::vector<T>& rOutput)
{
    std::ifstream fs(sPathToFile, std::ios::in | std::ios::binary);
    if (fs.is_open())
    {
        fs.seekg(0, std::ios::end);
        size_t nFileSize = fs.tellg();
        fs.seekg(0, std::ios::beg);

        rOutput.resize(nFileSize / sizeof(T));
        fs.read(reinterpret_cast<char*>(rOutput.data()), nFileSize);
        fs.close();
        return true;
    }
    else
    {
        return false;
    }
}

std::string getExeFolder()
{
    CHAR path[MAX_PATH];
    if (GetModuleFileNameA(NULL, path, MAX_PATH) == 0)
    {
        return std::string();
    }
    std::string sExePath = path;
    return sExePath.substr(0, sExePath.find_last_of("\\/")) + "/";
}

int main()
{
    const std::string sExeFolder = getExeFolder();
    std::vector<std::string> tasks = { "afiro", "OSA-60" };
    for (auto& task : tasks)
    {
        ClpSimplex initialModel;
        int status = initialModel.readMps((sExeFolder + "data/mps/" + task + ".mps").c_str());
        if (status)
        {
            return status;
        }

        AgmonMotzkin agmonMotzkin;
        if (agmonMotzkin.init(initialModel) == false)
        {
            return -1;
        }

        std::string sInitialPointPath(sExeFolder + "data/feasible_solutions/" + task +".bin");
        if (fs::exists(sInitialPointPath))
        {
            std::vector<double> initialPoint;
            if (loadBinFile<double>(sInitialPointPath, initialPoint) == false)
            {
                return -1;
            }
            if (agmonMotzkin.setCurrentPoint(initialPoint) == false)
            {
                return -1;
            }
        }
        {
            Timer t("Agmon-Motzkin algorithm time: ");
            while (agmonMotzkin.nextApproximation(1e-8));
        }
        double dObjectiveValue = agmonMotzkin.getCurrentObjectiveValue();
        printf("Agmon-Motzkin algorithm objective value = %f\n", dObjectiveValue);

        auto initialSolution = agmonMotzkin.getCurrentPoint();
        if (saveToBinFile<double>(sExeFolder + task + ".bin", initialSolution) == false)
        {
            std::cout << "Can not save initial solution to file\n";
        }

        std::ofstream resultFile(task + "_result.csv", std::ios::out | std::ios::trunc);
        if (resultFile.is_open() == false)
        {
            return -1;
        }
        resultFile << "Model,Solve Type,Presolve Type,Iterations,Time (s)\n";
        for (int iSolveType = 0; iSolveType < 6; ++iSolveType)
        {
            ClpSolve solverOptions;
            solverOptions.setSolveType(ClpSolve::SolveType(iSolveType));
            for (int iPresolveType = 0; iPresolveType < 4; ++iPresolveType)
            {
                solverOptions.setPresolveType(ClpSolve::PresolveType(iPresolveType));
                std::cout << '\n' << "usualModel:\n";
                auto fromBegin = initialModel;
                double fromBeginTime = 0;
                {
                    Timer t("Solving time: ");
                    fromBegin.initialSolve(solverOptions);
                    fromBeginTime = t.getSeconds();
                }
                int fromBeginIterations = fromBegin.getIterationCount();
                std::cout << "Iterations count: " << fromBeginIterations << '\n';

                std::cout << '\n' << "fromFeasible:\n";
                auto fromFeasible = initialModel;
                fromFeasible.setColSolution(initialSolution.data());
                double fromFeasibleTime = 0;
                {
                    Timer t("Solving time: ");
                    fromFeasible.initialSolve(solverOptions);
                    fromFeasibleTime = t.getSeconds();
                }
                int fromFeasibleIterations = fromFeasible.getIterationCount();
                std::cout << "Iterations count: " << fromFeasibleIterations << '\n';

                resultFile << "From begin," << iSolveType << ',' << iPresolveType << ','
                    << fromBeginIterations << ',' << fromBeginTime << '\n';
                resultFile << "From feasible," << iSolveType << ',' << iPresolveType << ','
                    << fromFeasibleIterations << ',' << fromFeasibleTime << '\n';
            }
        }
        resultFile.close();
    }
    return 0;
}
