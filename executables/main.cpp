#include <iostream>
#include <numeric>
#include <fstream>
#include <filesystem>

#include <AgmonMotzkin.h>
#include <File.h>
#include <Timer.h>

namespace fs = std::filesystem;

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
            constexpr double dEpsilon = 1e-8;
            constexpr size_t nMaxSteps = 1000000;
            size_t nStep = 0;
            while (++nStep < nMaxSteps)
            {
                auto [iRow, dStep] = agmonMotzkin.prepareStep();
                if ((iRow < 0) || (dStep < dEpsilon))
                {
                    break;
                }
                agmonMotzkin.updateCurrentPoint(iRow, dStep);
            }
        }
        double dObjectiveValue = agmonMotzkin.getCurrentObjectiveValue();
        std::cout << "Agmon-Motzkin algorithm objective value = " << dObjectiveValue << '\n';

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
