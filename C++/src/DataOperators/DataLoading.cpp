#include "DataLoading.h"
std::vector<int> Bins;
std::vector<std::string> Files;
std::vector<int> NumberOfStarsInFile;

std::vector<int> StarsLeftInFile;
std::vector<std::vector<int>> batchCounts;
std::vector<std::vector<int>> batchStarts;
std::vector<int> batchOffsets;

void GetAssignments(int id,std::string dataSource)
{
	//collect number of stars per file
	std::string starDirectoryFile = dataSource + "/directory.dat";
	std::vector<FileStarPairs> files;
	forLineVectorInFile(starDirectoryFile,' ',
			FileStarPairs f;
			f.FileName = FILE_LINE_VECTOR[1];
			f.NStars = std::stoi(FILE_LINE_VECTOR[0]);
			files.push_back(f);
	);
	
	
	std::string assignmentFile = "coreAssignments.dat";
	
	forLineVectorInFile(assignmentFile,',',
		
		int core = stoi(FILE_LINE_VECTOR[0]);
		if (core == id)
		{
			for (int i = 1; i < FILE_LINE_VECTOR.size(); i+=2)
			{
				std::string filename = FILE_LINE_VECTOR[i];
				Files.push_back(dataSource + filename);
				Bins.push_back(stoi(FILE_LINE_VECTOR[i+1]));
				
				bool fileInDirectory = false;
				for (int j = 0; j < files.size(); ++j)
				{
					if (filename == files[j].FileName)
					{
							int n = files[j].NStars;
							
							if (DataLoadCount > 0)
							{
								n = std::min(n,DataLoadCount);
							}
							NumberOfStarsInFile.push_back(n);
							StarsLeftInFile.push_back(n);
							fileInDirectory = true;
							break;
					}
				}
				if (!fileInDirectory)
				{
					ERROR(2,"A file was in the core allocation file, but was not included in the datasource directory");
				}
				
			}
		}
	);
}


void CalculateBatches(int id)
{
	int bRemaining = N_SGD_Batches;
	int NFiles = Files.size();
	batchCounts = std::vector<std::vector<int>>(N_SGD_Batches,std::vector<int>(NFiles,0));
	batchStarts = std::vector<std::vector<int>>(N_SGD_Batches,std::vector<int>(NFiles,0));
	batchOffsets = std::vector<int>(N_SGD_Batches,0);
	int batchOffset = 0;
	
	while (bRemaining > 1)
	{
		
		int index = N_SGD_Batches - bRemaining;
		int starsRemaining = std::accumulate(StarsLeftInFile.begin(),StarsLeftInFile.end(),0);
		int starsPerBatch = starsRemaining / bRemaining;
		
		int sum = 0;
		for (int i = 0; i < NFiles; ++i)
		{
			int contribution = starsPerBatch*(float)StarsLeftInFile[i]/starsRemaining;
			batchCounts[index][i] = contribution;
			StarsLeftInFile[i] -= contribution;
			batchStarts[index][i] = sum + batchOffset;
			sum += contribution;
		}

		batchOffsets[index] = batchOffset;
		batchOffset += sum;
		--bRemaining; 
	}
	
	int finalBatchSize = 0;
	for (int i = 0; i < NFiles; ++i)
	{
		batchStarts[N_SGD_Batches - 1][i] = batchOffset + finalBatchSize;
		finalBatchSize += StarsLeftInFile[i];
		batchCounts[N_SGD_Batches - 1][i] = StarsLeftInFile[i];
		
	}
	batchOffsets[N_SGD_Batches - 1] = batchOffset;
	
	if (id == RootID)
	{
		std::cout << "BATCH ALLOCATIONS!" << std::endl;
		for (int i = 0; i < NFiles; ++i)
		{
			for (int j = 0; j < N_SGD_Batches; ++j)
			{
				std::cout <<std::setw(7) <<batchCounts[j][i]; 
			}
			std::cout << "\n";
		}
	}
	
	
	

}

std::vector<int>  LoadData(const int ProcessRank, const int JobSize, std::vector<Star> & Data, int & TotalStars,const std::string dataSource)
{
	if (ProcessRank == RootID)
	{
		GlobalLog(1,
			std::cout << "Initialising starAllocation script...\n";
			std::string command = "python3 starAllocation.py " + dataSource + " " + std::to_string(JobSize);
			std::cout << command << std::endl;
			system(command.c_str() );
			std::cout << "Data allocation complete, beginning readin...." <<std::endl;
		);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	auto start = std::chrono::system_clock::now();
	GetAssignments(ProcessRank,dataSource);
	CalculateBatches(ProcessRank);
	int nStarsAssigned = std::accumulate(NumberOfStarsInFile.begin(),NumberOfStarsInFile.end(),0);
	
	Data.resize(nStarsAssigned);
	
	
	bool isReporter = (ProcessRank == JobSize - 1);
	int meaningfullyLargeNumber = 1e8;
	int readIn = 0;
	int lastCheckPoint = 0;
	
	int fileBatchShift = 0;
	for (int i = 0; i < Files.size(); ++i)
	{
		std::string file = Files[i];
		int gBin = Bins[i];

		//use a fancy macro (FileHandler.h)
		int idx = 0;
		int batch = 0;
		while (batchCounts[batch][i] == 0)
		{
			++batch;
		}
		int batchOffset = batchStarts[0][i];
		int starsLoaded = 0;
		forLineVectorInFile(file,',',
			
			int loc = batchStarts[batch][i] + idx;
			Data[loc] = Star(FILE_LINE_VECTOR,gBin);
			std::cout << "Loaded a star with " << Data[loc].nMeasure << "  " << Data[loc].nVisit << "  " << Data[loc].TimeSeries.size() << "  into batch  " << batch << std::endl;
			 
			++starsLoaded;
			++idx;
			
			if (idx == batchCounts[batch][i])
			{
				++batch;
				idx = 0;
			}
			
			if (DataLoadCount > 0 && starsLoaded >= DataLoadCount)
			{
				break;
			}
		);
	}

	GlobalLog(1,
		auto end = std::chrono::system_clock::now();
		std::string duration = formatDuration(start,end);
		std::cout << "\tProcess " << ProcessRank << " has loaded in " << Data.size() << " datapoints in " << duration << std::endl; 
	);
	
	int n = Data.size();

	int MaxStarsInCore = 0;
	MPI_Reduce(&n,&TotalStars,1,MPI_INT,MPI_SUM,RootID,MPI_COMM_WORLD);
	MPI_Reduce(&n, &MaxStarsInCore, 1,MPI_INT, MPI_MAX, RootID,MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (ProcessRank==RootID)
	{
		GlobalLog(0,
			std::cout << TotalStars << " stars have been loaded into memory (max stars in core: " << MaxStarsInCore << ")" << std::endl;
		);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	return batchOffsets;
}

