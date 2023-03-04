#include "cpuTime.H" //OF class for timing, source file located in /src/OSspecific/POSIX/cpuTime
cpuTime timer;

List<scalar> nonReactingTime; //variables
List<scalar> reactingTime;
List<scalar> matrixTime;
nonReactingTime.setSize(Pstream::nProcs(),-VGREAT);
reactingTime.setSize(Pstream::nProcs(),-VGREAT);
matrixTime.setSize(Pstream::nProcs(),-VGREAT);
nonReactingTime[Pstream::myProcNo()]=0.0;
reactingTime[Pstream::myProcNo()]=0.0;
matrixTime[Pstream::myProcNo()]=0.0;

enum timerType {nonReacting, reacting, matrix}; //timing strategy
timerType currentTimerType=nonReacting;

std::fstream cpuTiming_file; //output
