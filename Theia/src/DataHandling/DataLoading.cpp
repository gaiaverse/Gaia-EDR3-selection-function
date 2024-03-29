#include "DataLoading.h"
std::vector<int> Bins;

std::vector<int> NumberOfStarsInFile;
std::vector<int> StarsLeftInFile;
std::vector<std::vector<int>> batchCounts;


std::vector<File> GetAssignments(int id,std::string dataSource)
{
	//collect number of stars per file
	std::string starDirectoryFile = dataSource + "/directory.dat";
	std::vector<FileStarPairs> files;
	forLineVectorIn(starDirectoryFile,' ',
			FileStarPairs f;
			f.FileName = FILE_LINE_VECTOR[1];
			f.NStars = std::stoi(FILE_LINE_VECTOR[0]);
			files.push_back(f);
	);
	std::vector<File> Files;
	
	std::string assignmentFile = "../../ModelInputs/coreAssignments.dat";
	
	forLineVectorIn(assignmentFile,',',
		
		int core = stoi(FILE_LINE_VECTOR[0]);
		if (core == id)
		{
			for (int i = 1; i < FILE_LINE_VECTOR.size(); i+=2)
			{
				std::string filename = FILE_LINE_VECTOR[i];
				
				File f;
				f.Name = dataSource + filename;
				f.Bin = stoi(FILE_LINE_VECTOR[i+1]) - magOffset;
				
				Files.push_back(f);
				
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
	return Files;
}


void CalculateBatches(std::vector<File> & Files, int batches)
{
	
	
	int bRemaining = batches;
	int NFiles = Files.size();
	batchCounts = std::vector<std::vector<int>>(batches,std::vector<int>(NFiles,0));
	int batchOffset = 0;
	
	while (bRemaining > 1)
	{
		
		int index = batches - bRemaining;
		int starsRemaining = std::accumulate(StarsLeftInFile.begin(),StarsLeftInFile.end(),0);
		int starsPerBatch = starsRemaining / bRemaining;
		
		int sum = 0;
		for (int i = 0; i < NFiles; ++i)
		{
			int contribution = starsPerBatch*(float)StarsLeftInFile[i]/starsRemaining;
			batchCounts[index][i] = contribution;
			StarsLeftInFile[i] -= contribution;
			sum += contribution;
		}

		--bRemaining; 
	}
	
	int finalBatchSize = 0;
	for (int i = 0; i < NFiles; ++i)
	{
		finalBatchSize += StarsLeftInFile[i];
		batchCounts[batches - 1][i] = StarsLeftInFile[i];	
	}
}

void  LoadData(const int ProcessRank, const int JobSize, std::vector<std::vector<Star>> & Data, int & TotalStars,const std::string dataSource, int batches)
{
	if (ProcessRank == RootID)
	{	
		std::cout << "Initialising starAllocation script...\n";
		std::string command = "python3 src/DataHandling/starAllocation.py " + dataSource + " " + std::to_string(JobSize);
		system(command.c_str() );
		std::cout << "Data allocation complete, beginning readin...." <<std::endl;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	auto start = std::chrono::system_clock::now();
	std::vector<File> Files = GetAssignments(ProcessRank,dataSource);
	CalculateBatches(Files,batches);

	//resize datafile
	Data.resize(batches);

	int allStarsLoaded = 0;
	for (int i = 0; i < Files.size(); ++i)
	{
		std::string file = Files[i].Name;
		int gBin = Files[i].Bin;

		//use a fancy macro (FileHandler.h)
		
		int batch = 0;
		int idx =0;
		
		int starsLoaded = 0;
		forLineVectorIn(file,',',

			while (batchCounts[batch][i] == 0)
			{
				++batch;
			}

			Data[batch].push_back(Star(FILE_LINE_VECTOR,gBin));

			++starsLoaded;
			++allStarsLoaded;
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
	auto end = std::chrono::system_clock::now();
	std::string duration = JSL::FormatTimeDuration(start,end);
	std::cout << "\tProcess " << ProcessRank << " has loaded in " << allStarsLoaded << " datapoints in " << duration << std::endl; 

	
	
	int MaxStarsInCore = 0;
	MPI_Reduce(&allStarsLoaded,&TotalStars,1,MPI_INT,MPI_SUM,RootID,MPI_COMM_WORLD);
	MPI_Reduce(&allStarsLoaded, &MaxStarsInCore, 1,MPI_INT, MPI_MAX, RootID,MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (ProcessRank==RootID)
	{
		std::cout << TotalStars << " stars have been loaded into memory (max stars in core: " << MaxStarsInCore << ")" << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

