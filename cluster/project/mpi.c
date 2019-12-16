#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <sys/time.h>
#define LINE 100000

void readfile(double num[LINE]);
void printfile(double num[LINE]);
void copyData(double num[LINE], double num1[LINE]);
int merge(double *ina, int lena, double *inb, int lenb, double *out);
int domerge_sort(double *a, int start, int end, double *b);
int merge_sort(int n, double *a);
void MPI_Pairwise_Exchange(int localn, double *locala, int sendrank, int recvrank, MPI_Comm comm);
void mpi(double num[LINE], double num1[LINE]);

int main(int argc, char **argv) {
	int rank;
	double num[LINE];
	double num1[LINE];
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank == 0){
		printf("Getting data...\n");
		readfile(num);
		printf("Sorting data...\n\n");
	}
	mpi(num, num1);
	if(rank == 0){
		printfile(num1);
		printf("\nParallel odd even sort in mpi sucessfully.\n");
	}
    MPI_Finalize();
    return 0;
}

void readfile(double num[LINE]){
	double temp;
	int i;
	FILE *fp;
	fp = fopen("number.txt", "r");
	i = 0;
	if(fp == NULL){
		printf("Error loading file!!\n");
		exit(1);
	}else{
		while(!feof(fp)){
			fscanf(fp, "%lf", &temp);
			num[i] = temp;
			i++;
		}
	}
	fclose(fp);
}

void printfile(double num[LINE]){
	int i;
	FILE *fp = fopen("update.txt", "w");
	for (i = 0; i < LINE; i++)
		fprintf(fp, "%lf   ", num[i]);
	fclose(fp);
}

void copyData(double num[LINE], double num1[LINE]){
	int i;
	for(i = 0; i < LINE; i++)
		num1[i] = num[i];
}

int merge(double *ina, int lena, double *inb, int lenb, double *out){
    int i,j;
    int outcount=0;

    for (i=0,j=0; i<lena; i++) {
        while ((inb[j] < ina[i]) && j < lenb) {
            out[outcount++] = inb[j++];
        }
        out[outcount++] = ina[i];
    }
    while (j<lenb)
        out[outcount++] = inb[j++];

    return 0;
}

int domerge_sort(double *a, int start, int end, double *b){
	int i;
    if ((end - start) <= 1) return 0;

    int mid = (end+start)/2;
    domerge_sort(a, start, mid, b);
    domerge_sort(a, mid,   end, b);
    merge(&(a[start]), mid-start, &(a[mid]), end-mid, &(b[start]));
    for (i=start; i<end; i++)
        a[i] = b[i];

    return 0;
}

int merge_sort(int n, double *a){
    double b[n];
    domerge_sort(a, 0, n, b);
    return 0;
}

void MPI_Pairwise_Exchange(int localn, double *locala, int sendrank, int recvrank, MPI_Comm comm){
    int rank, i;
    double remote[localn];
    double all[2*localn];
    const int mergetag = 1;
    const int sortedtag = 2;

    MPI_Comm_rank(comm, &rank);
    if (rank == sendrank) {
        MPI_Send(locala, localn, MPI_DOUBLE, recvrank, mergetag, MPI_COMM_WORLD);
        MPI_Recv(locala, localn, MPI_DOUBLE, recvrank, sortedtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(remote, localn, MPI_DOUBLE, sendrank, mergetag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        merge(locala, localn, remote, localn, all);

        int theirstart = 0, mystart = localn;
        if (sendrank > rank) {
            theirstart = localn;
            mystart = 0;
        }
        MPI_Send(&(all[theirstart]), localn, MPI_DOUBLE, sendrank, sortedtag, MPI_COMM_WORLD);
        for (i=mystart; i<mystart+localn; i++)
            locala[i-mystart] = all[i];
    }
}

int MPI_OddEven_Sort(int n, double *a, int root, MPI_Comm comm){
    int rank, size, i;
    double *local_a;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    local_a = (double *) calloc(n / size, sizeof(double));
    MPI_Scatter(a, n / size, MPI_DOUBLE, local_a, n / size, MPI_DOUBLE, root, comm);
    merge_sort(n / size, local_a);
    for(i = 1; i <= size; i++){
        if((i + rank) % 2 == 0){  // means i and rank have same nature
            if(rank < size - 1){
                MPI_Pairwise_Exchange(n / size, local_a, rank, rank + 1, comm);
            }
        }else if(rank > 0){
            MPI_Pairwise_Exchange(n / size, local_a, rank - 1, rank, comm);
        }
    }
    MPI_Gather(local_a, n / size, MPI_DOUBLE, a, n / size, MPI_DOUBLE, root, comm);
    return MPI_SUCCESS;
}

void mpi(double num[LINE], double num1[LINE]){
	int line, rank, i;
	struct timeval tv;
	struct timezone tz;
	double start, end, time, time1, time2, average;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	start = 0;
	end = 0;
	time = 0;
	time1 = 0;
	time2 = 0;
	line = 10000;
	average = 0;
	if(rank == 0){
		printf("Time execution for parallel odd even sort using mpi\n");
		printf("================================================================================\n");
		printf("   Number of data        1st time       2nd time       3rd time       average   \n");
		printf("================================================================================\n");
	}
	while (line <= LINE){
		for (i = 0; i < 3; i++){
			if(rank == 0){
				copyData(num, num1);
				gettimeofday(&tv, &tz);
				start = (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
			}
			MPI_OddEven_Sort(line, num1, 0, MPI_COMM_WORLD);
			if(rank == 0){
				gettimeofday(&tv, &tz);
				end = (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
		
				if (i == 0)
					time = end - start;
				else if (i == 1)
					time1 = end - start;
				else if (i == 2)
					time2 = end - start;\
			}
		}
		if(rank == 0){
			average = (time + time1 + time2) / 3;
			printf("      %i              %fs      %fs      %fs      %fs\n", line, time, time1, time2, average);
		}
		line += 10000;
	}
}
