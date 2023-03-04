for (int i=0; i<Pstream::nProcs(); i++)
{
	reduce(nonReactingTime[i],maxOp<scalar>());
	reduce(reactingTime[i],maxOp<scalar>());
	reduce(matrixTime[i],maxOp<scalar>());
}

if(Pstream::master())
{
	cpuTiming_file.open("cpuTiming.txt", std::fstream::out | std::fstream::trunc);
	cpuTiming_file << "ProcessorID non_reacting_time reacting_time communication_and_matrix_operation_time" << std::endl;
	for (int i=0; i< Pstream::nProcs(); i++)
	{
		cpuTiming_file << i << " " << nonReactingTime[i] << " " << reactingTime[i] << " " << matrixTime[i] << std::endl;
	}
	cpuTiming_file.close();
}
