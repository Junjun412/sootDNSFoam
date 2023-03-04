switch (currentTimerType)
{
	case nonReacting:
		nonReactingTime[Pstream::myProcNo()] += timer.cpuTimeIncrement(); //time increment since last mark is added into a certain type
		break;
	case reacting:
		reactingTime[Pstream::myProcNo()] += timer.cpuTimeIncrement();
		break;
	case matrix:
		matrixTime[Pstream::myProcNo()] += timer.cpuTimeIncrement();
		break;
}
