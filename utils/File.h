#ifndef MV_File_H_

#include <Windows.h>

#include <fstream>
#include <vector>

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

#endif // MV_File_H_