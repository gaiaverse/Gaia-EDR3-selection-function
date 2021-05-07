#include "DataLoading.h"
std::vector<int> Bins;
std::vector<std::string> Files;
std::vector<int> NumberOfStarsInFile;

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
							NumberOfStarsInFile[j] = files[j].NStars;
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

void LoadData(const int ProcessRank, const int JobSize, std::vector<Star> * Data, int & TotalStars,std::string dataSource)
{
	if (ProcessRank == RootID)
	{
		GlobalLog(1,
			std::cout << "Initialising starAllocation script...\n";
			std::string command = "python starAllocation.py " + dataSource + " " + std::to_string(JobSize);
			system(command.c_str() );
			std::cout << "Data allocation complete, beginning readin...." <<std::endl;
		);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	auto start = std::chrono::system_clock::now();
	GetAssignments(ProcessRank,dataSource);
	bool isReporter = (ProcessRank == JobSize - 1);
	int meaningfullyLargeNumber = 1e8;
	int readIn = 0;
	int lastCheckPoint = 0;
	
	for (int i = 0; i < Files.size(); ++i)
	{
		std::string file = Files[i];
		int gBin = Bins[i];
		//use a fancy macro (FileHandler.h)
		forLineVectorInFile(file,',',
			Star s = Star(FILE_LINE_VECTOR,gBin);
			Data->push_back(s);
		);
	}
	
	
	
	GlobalLog(1,
		auto end = std::chrono::system_clock::now();
		std::string duration = formatDuration(start,end);
		std::cout << "\tProcess " << ProcessRank << " has loaded in " << Data->size() << " datapoints in " << duration << std::endl; 
	);
	
	int n = Data->size();

	int MaxStarsInCore = 0;
	MPI_Reduce(&n,&TotalStars,1,MPI_INT,MPI_SUM,RootID,MPI_COMM_WORLD);
	MPI_Reduce(&n, &MaxStarsInCore, 1,MPI_INT, MPI_MAX, RootID,MPI_COMM_WORLD);
	
	if (ProcessRank==RootID)
	{
		GlobalLog(0,
			std::cout << TotalStars << " stars have been loaded into memory (max stars in core: " << MaxStarsInCore << ")" << std::endl;
		);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

