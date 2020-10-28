#include <mpi.h>

#include <AgmonMotzkin.h>
#include <File.h>

int main(int argc, char* argv[])
{
    int iStatus = MPI_Init(&argc, &argv);
    if (iStatus != MPI_SUCCESS)
    {
        return MPI_Abort(MPI_COMM_WORLD, iStatus);
    }
    int iNumProcesses = -1;
    (void)MPI_Comm_size(MPI_COMM_WORLD, &iNumProcesses);
    int iProcId = -1;
    (void)MPI_Comm_rank(MPI_COMM_WORLD, &iProcId);

    const std::string sExeFolder = getExeFolder();
    std::vector<std::string> tasks = { "afiro" };
    for (auto& task : tasks)
    {
        ClpSimplex initialModel;
        int iStatus = initialModel.readMps((sExeFolder + "data/mps/" + task + ".mps").c_str());
        if (iStatus)
        {
            return iStatus;
        }

        AgmonMotzkin agmonMotzkin;
        if (agmonMotzkin.init(initialModel) == false)
        {
            return -1;
        }

        double dStartWtime = 0;
        if (iProcId == 0)
        {
            dStartWtime = MPI_Wtime();
        }

        constexpr size_t nIntSize = sizeof(int);
        constexpr size_t nDoubleSize = sizeof(double);
        constexpr size_t nTotalSize = nIntSize + nDoubleSize;
        std::vector<char> buffer(nTotalSize);
        int* pRow = reinterpret_cast<int*>(buffer.data());
        double* pStep = reinterpret_cast<double*>(buffer.data() + nIntSize);
        constexpr int iTag = 0;

        constexpr double dEpsilon = 1e-8;
        constexpr size_t nMaxSteps = 1000000;
        size_t nStep = 0;
        while (++nStep < nMaxSteps)
        {
            auto [iRow, dStep] = agmonMotzkin.prepareStep(iProcId, iNumProcesses);
            
            if (iProcId == 0)
            {
                for (int iId = 1; iId < iNumProcesses; ++iId)
                {
                    MPI_Status status;
                    (void)MPI_Recv(buffer.data(), nTotalSize, MPI_CHAR, iId, iTag, MPI_COMM_WORLD, &status);
                    if (dStep < *pStep)
                    {
                        dStep = *pStep;
                        iRow = *pRow;
                    }
                }
                *pStep = dStep;
                *pRow = iRow;
            }
            else
            {
                *pRow = iRow;
                *pStep = dStep;
                (void)MPI_Send(buffer.data(), nTotalSize, MPI_CHAR, 0, iTag, MPI_COMM_WORLD);
            }
            (void)MPI_Bcast(buffer.data(), nTotalSize, MPI_CHAR, 0, MPI_COMM_WORLD);
            if ((*pRow < 0) || (*pStep < dEpsilon))
            {
                break;
            }
            agmonMotzkin.updateCurrentPoint(*pRow, *pStep);
        }

        if (iProcId == 0)
        {
            double dObjectiveValue = agmonMotzkin.getCurrentObjectiveValue();
            printf("Agmon-Motzkin algorithm objective value = %f\n", dObjectiveValue);

            double dEndWtime = MPI_Wtime();
            printf("Calculation time: %f s\n", (dEndWtime - dStartWtime));

            auto initialSolution = agmonMotzkin.getCurrentPoint();
            if (saveToBinFile<double>(sExeFolder + task + ".bin", initialSolution) == false)
            {
                printf("Can not save initial solution to file\n");
            }
        }
    }

    MPI_Finalize();
    return 0;
}