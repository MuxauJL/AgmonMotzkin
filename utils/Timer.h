#ifndef MV_Timer_H_
#include <chrono>
#include <iostream>

template<typename clock = std::chrono::high_resolution_clock, typename time_point = clock::time_point>
class Timer
{
public:
    Timer(const std::string& sMessage = "")
        : sMessage_(sMessage), startTime_(clock::now())
    {}

    ~Timer()
    {
        std::cout << sMessage_
            << std::chrono::duration_cast<std::chrono::microseconds>(clock::now() - startTime_).count()
            << " microseconds" << std::endl;
    }

    inline double getSeconds()
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(clock::now() - startTime_).count() / 1e6;
    }

    inline double getMilliseconds()
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(clock::now() - startTime_).count() / 1e3;
    }

    inline double getMicroseconds()
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(clock::now() - startTime_).count();
    }
private:
    std::string sMessage_;
    time_point startTime_;
};

#endif // MV_Timer_H_
