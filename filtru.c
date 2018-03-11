#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <math.h>
#include<mpi.h>

#define L 100

MPI_Status status;
MPI_Request request;
	
void readLine(char* line, int lineNumber, char filename[]) {

	FILE *file = fopen(filename, "r");
	int count = 0;
	if ( file != NULL ) {

	    while (fgets(line, L, file) != NULL) {

	        if (count == lineNumber) {
	        	return;
	        }
	        else {
	            count++;
	        }
	    }
	    fclose(file);
	}
}

void makeNode(int* children, int parent, int rank, char* line) {
	int child, i, nada;
	char *p, piw[2] = " ";
	p = strtok(line, piw);

	i = 0;	
	while (p != NULL) {

		p = strtok(NULL, piw);
		if (p==NULL)
			break;

		child = atoi(p);

		if (child != parent) {
			children[i++] = child;
			MPI_Send(&rank, 1, MPI_INT, child, 0, MPI_COMM_WORLD);
		}
	}
}

void selectPicture(char* line, char* image, char* out_image, int* filter) {

	char *p, piw[2] = " ";
	p = strtok(line, piw);
	if (strcmp("sobel", p) == 0) *filter = 1;
	if (strcmp("mean_removal", p) == 0) *filter = 2;
	p = strtok(NULL, piw);
	strcpy(image, p);

	p = strtok(NULL, piw);
	strcpy(out_image, p);
	out_image[strlen(out_image)-1] = '\0';
}

int* uploadPicture(char* filename, int* height, int* width, char* line1, char* line2, char* line4) {

	FILE *file = fopen(filename, "r");
	char *line, piw[2] = " ", *p;
	int *matrix = NULL, i = 0, k;

	if ( file != NULL ) {
	
		fgets(line1, L, file);
		fgets(line2, L, file);
		fgets(line, L, file);
				
		p = strtok(line, piw);
		*width = atoi(p) + 2;

		p = strtok(NULL, piw);
		*height = atoi(p);

		fgets(line4, L, file);

		matrix = (int *) malloc((*height) * (*width) * sizeof(int));

	    while (fgets(line, L, file) != NULL) {
			
			if (i%(*width) == 0)
				matrix[i++] = 0;

			matrix[i++] = atoi(line);

	    	if ((i+1)%(*width)==0)
				matrix[i++] = 0;
	    }

	}
    fclose(file);

	return matrix;
}

void sendBlocks(int* children, int no_of_children, int height, int width, int* matrix, int filter, int* Umargin, int* Dmargin, int k) {

	int no_of_lines = height/no_of_children;
	int rest = height%no_of_children;

	int i = 0;
	int umargin[width], dmargin[width];

	while (i < no_of_children) {
		int ii;

		if (i == no_of_children -1) 
			no_of_lines += rest;

		MPI_Send(&width, 1, MPI_INT, children[i], filter, MPI_COMM_WORLD);
		MPI_Send(&no_of_lines, 1, MPI_INT, children[i], filter, MPI_COMM_WORLD);

		if (i == 0) {
			for (ii = 0; ii < width; ii++)
				umargin[ii] = Umargin[ii];
		}
		else {
			for (ii = 0; ii < width; ii++)
				umargin[ii] = matrix[(i*no_of_lines-1)*width + ii];
		} 

		MPI_Send(umargin, width, MPI_INT, children[i], filter, MPI_COMM_WORLD);

		if (i == no_of_children - 1) {
			for (ii = 0; ii < width; ii++)
				dmargin[ii] = Dmargin[ii];
		}
		else {
			for (ii = 0; ii < width; ii++) 
				dmargin[ii] = matrix[(i+1)*no_of_lines*width + ii];
		}
		MPI_Send(dmargin, width, MPI_INT, children[i], filter, MPI_COMM_WORLD);
			
		
		MPI_Send(matrix + i*width*no_of_lines, no_of_lines * width, MPI_INT, children[i], filter, MPI_COMM_WORLD);

		i++;
	}
} 

int* recvMargins(int rank, int k, int parent, int *filter, int *height, int *width) {
	
	int aux, i, *AUX;
	int *Margins = NULL;

	MPI_Recv(&aux, 1, MPI_INT, parent, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	*width = aux;

	*filter = status.MPI_TAG;
	MPI_Recv(&aux, 1, MPI_INT, parent, *filter, MPI_COMM_WORLD, &status);
	*height = aux;
	
	Margins = (int *) malloc(2*(*width)*sizeof(int));
	AUX = (int *) malloc((*width)*sizeof(int));

	MPI_Recv(AUX, *width, MPI_INT, parent, *filter, MPI_COMM_WORLD, &status);

	for (i = 0; i < (*width); i++)
		Margins[i] = AUX[i];

	MPI_Recv(AUX, *width, MPI_INT, parent, *filter, MPI_COMM_WORLD, &status);
	for (i = 0; i < (*width); i++) 
		Margins[(*width)+i] = AUX[i];

	free(AUX);
	return Margins;
}

int* recvBlock(int parent, int filter, int height, int width) {

	int *matrix = NULL, i = 0;
	matrix = (int *) malloc(height*width*sizeof(int));
	
	MPI_Recv(matrix, height * width, MPI_INT, parent, filter, MPI_COMM_WORLD, &status);

	return matrix;
}

int* applyFilter(int* matrix, int height, int width, int* Umargin, int* Dmargin, int** sepia, int factor, int displ) {

	int *filtered_img = NULL, k,i;
	filtered_img = (int *) malloc (width * height * sizeof(int));

	for (k = 0; k < height; k++) {
		for (i = 1; i < width-1; i++) {
			int ii, kk;
			float val = 0;

			for (kk = 0; kk < 3; kk++) {
				for (ii = 0; ii < 3; ii++) {

					if (k == 0 && kk == 0) {
						val += sepia[kk][ii] * Umargin[i+ii-1];
						continue;
					}

					if (k == height-1 && kk == 2) {
						val += sepia[kk][ii] * Dmargin[i+ii-1];
						continue;
					}

					val += sepia[kk][ii] * matrix[(k+kk-1)*width+i+ii-1];
				}
			}

			val /= factor;
			val += displ;

			if (val - floor(val) < 0.5)
				filtered_img[k*width+i] = floor(val);
			else
				filtered_img[k*width+i] = floor(val+1);

			if (val < 0) 
				filtered_img[k*width+i] = 0;
			if (val > 255) 
				filtered_img[k*width+i] = 255;

		}
	}

	return filtered_img;
}

int* recvImage(int* children, int no_of_children, int height, int width, int rank) {
	int no_of_lines = height/no_of_children;
	int rest = height%no_of_children;
	int *img = NULL;
	img = (int *) malloc (height * width * sizeof(int));

	int i = 0;

	while (i < no_of_children) {
		int *AUX = NULL, ii;
		if (i == no_of_children -1) 
			no_of_lines += rest;
		AUX = (int *) malloc(no_of_lines * width * sizeof(int));
		MPI_Recv(AUX, no_of_lines*width, MPI_INT, children[i], 0, MPI_COMM_WORLD, &status);

		for (ii = 0; ii < no_of_lines * width; ii++)
			img[i * no_of_lines * width + ii] = AUX[ii];
		free(AUX);
		i++;
	}

	return img;
}

void wirteM2file(char* out_image, char* line1, char* line2, char* line4, int height, int width, int* matrix) {
	FILE *file = fopen(out_image, "w");
	char aux[10];
	int i, k, w = width, h = height;

	if ( file != NULL ) {
		fputs(line1, file);
		fputs(line2, file);

		sprintf(aux, "%d ", w-2);
		fputs(aux, file);
		
		sprintf(aux, "%d\n", h);
		fputs(aux, file);
		
		fputs(line4, file);

		for (k = 0; k < height; k++) {
 			for (i = 1; i < width-1; i++) {
				sprintf(aux, "%d\n", matrix[k*width+i]);
				fputs(aux, file);
			}
		}
		fclose(file);
	}
}

void recvStats(int* children, int no_of_children, int* statistics, int dim) {
	int i;
	
	for (i = 0; i < no_of_children; i++) {
		int *s ,ii;
		s = (int *) malloc(dim * sizeof(int));
		MPI_Recv(s, dim, MPI_INT, children[i], 0, MPI_COMM_WORLD, &status);
		for (ii = 0; ii < dim; ii++) 
			statistics[ii] |= s[ii];
		free(s);
	}

}

void writeS2file(char* filename, int* statistics, int dim) {
	FILE *file = fopen(filename, "w");
	char aux[10];
	int i;

	if (file != NULL) {

		for (i = 0; i < dim; i++) {
			sprintf(aux, "%d: %d\n", i, statistics[i]);
			fputs(aux, file);
		}
		fclose(file);
	}
}

int main(int argc, char * argv[]) {

	int rank, nProcesses;
	int parent = -1, nada, child, no_of_children;
	int i, k;
	int *matrix = NULL, *Dmargin = NULL, *Umargin = NULL;
	int width, height;
	char line[L], *p, piw[2] = " ";
	char line1[L], line2[L], line4[L];
	int *children = NULL, *statistics = NULL;
	int no_of_images, filter;
	char image[L], out_image[L];

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

	statistics = (int *) malloc(nProcesses * sizeof(int));
	for (i = 0; i < nProcesses; i++) {
		statistics[i] = 0;
	}

	if (rank == 0) {
		readLine(line, rank, argv[1]);
		
		k=0;		
		for (i=0; line[i+k]; line[i+k]==' ' ? i++ : k++);
		children = (int *) malloc(i * sizeof(int));
		no_of_children = i;
		makeNode(children, parent, rank, line);

		readLine(line, rank, argv[2]);
		no_of_images = atoi(line);

		for (i = 0; i < no_of_children; i++)
		 	MPI_Send(&no_of_images, 1, MPI_INT, children[i], 0, MPI_COMM_WORLD);

		k=0;
		while (no_of_images > k) {
			readLine(line, k+1, argv[2]);
			k++;
			selectPicture(line, image, out_image, &filter);
			matrix = uploadPicture(image, &height, &width, line1, line2, line4);
			{
				int zeros[width];
				for (i = 0; i < width; i++)
					zeros[i] = 0;

				sendBlocks(children, no_of_children, height, width, matrix, filter, zeros, zeros, k);
			}
		
			free(matrix);

			matrix = recvImage(children, no_of_children, height, width, rank);

			wirteM2file(out_image, line1, line2, line4, height, width, matrix);

		}
		recvStats(children, no_of_children, statistics, nProcesses);
		writeS2file(argv[3], statistics, nProcesses);

	}
	else {
		int stats = 0;

		MPI_Recv(&parent, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&no_of_images, 1, MPI_INT, parent, 0, MPI_COMM_WORLD, &status);

		readLine(line, rank, argv[1]);
		
		k = 0;
		for (i=0; line[i+k]; line[i+k]==' ' ? i++ : k++);
		children = (int *) malloc(--i * sizeof(int));
		no_of_children = i;
		
		makeNode(children, parent, rank, line);
		for (i = 0; i < no_of_children; i++)
		 	MPI_Send(&no_of_images, 1, MPI_INT, children[i], 0, MPI_COMM_WORLD);

		k = 0;
		

		while (no_of_images > k) {
			int *aux, *filtered_img;
			k++;

			aux = recvMargins(rank, k, parent, &filter, &height, &width);
			Umargin = (int *) malloc(width*sizeof(int));
			Dmargin = (int *) malloc(width*sizeof(int));

			for (i = 0; i < width; i++) {
				Umargin[i] = aux[i];
				Dmargin[i] = aux[i+width];
			}

			matrix = recvBlock(parent, filter, height, width);

			if (no_of_children != 0) {
				sendBlocks(children, no_of_children, height, width, matrix, filter, Umargin, Dmargin, k);
				// asteapta matricea procesata de la frunze
				filtered_img = recvImage(children, no_of_children, height, width, rank);
			}
			else {					// FILTRU
				int **sepia, factor, displ;
				stats+=height;
				sepia = (int **) malloc(3* sizeof(int*));

				for (i = 0; i < 3; i++)
					sepia[i] = (int *) malloc(3 * sizeof(int));

				if(filter == 1) {
					sepia[0][0] = 1; sepia[0][1] = 0; sepia[0][2] = -1;
					sepia[1][0] = 2; sepia[1][1] = 0; sepia[1][2] = -2;
					sepia[2][0] = 1; sepia[2][1] = 0; sepia[2][2] = -1;

					factor = 1;
					displ = 127;	
				} else {
					sepia[0][0] = -1; sepia[0][1] = -1; sepia[0][2] = -1;
					sepia[1][0] = -1; sepia[1][1] = 9; sepia[1][2] = -1;
					sepia[2][0] = -1; sepia[2][1] = -1; sepia[2][2] = -1;

					factor = 1;
					displ = 0;
				}

				filtered_img = applyFilter(matrix, height, width, Umargin, Dmargin, sepia, factor, displ);	
			}
			MPI_Send(filtered_img, height*width, MPI_INT, parent, 0, MPI_COMM_WORLD);
		}

		if (no_of_children != 0)
			recvStats(children, no_of_children, statistics, nProcesses);

		statistics[rank] = stats;
		MPI_Send(statistics, nProcesses, MPI_INT, parent, 0, MPI_COMM_WORLD);
	}

	free(children);
	MPI_Finalize();
	return 0;
}