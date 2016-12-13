/*
Fractal code for CS 4380 / CS 5351

Copyright (c) 2016, Texas State University. All rights reserved.

Redistribution in source or binary form, with or without modification,
is not permitted. Use in source and binary forms, with or without
modification, is only permitted for academic use in CS 4380 or CS 5351
at Texas State University.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Author: Martin Burtscher
editor: Diego Loya
*/

#include <iostream>
#include <mpi.h>
#include <math.h>
#include <cstdlib>
#include <sys/time.h>
#include "cs43805351.h"

static const double Delta = 0.005491;
static const double xMid = 0.745796;
static const double yMid = 0.105089;

int main(int argc, char *argv[])
{
  int my_rank;  //process number
  int comm_sz;  //number of processors
  int my_start; 
  int my_end;
  MPI_Init (NULL,NULL);  //all MPI function calls after this	
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  if (my_rank==0)
 	 printf("Fractal v1.5 [MPI]\n");




  // check command line
  if (argc != 3) {fprintf(stderr, "usage: %s frame_width num_frames\n", argv[0]); exit(-1);}
  int width = atoi(argv[1]);
  if (width < 10) {fprintf(stderr, "error: frame_width must be at least 10\n"); exit(-1);}
  int frames = atoi(argv[2]);
  if (frames < 1) {fprintf(stderr, "error: num_frames must be at least 1\n"); exit(-1);}
  if (my_rank==0)
  	printf("computing %d frames of %d by %d fractal\n", frames, width, width);

  // allocate picture array
  int array_size=frames/comm_sz;
    
  // allocation of picture array.  Process 0 will have a full size array, while the rest will only have a portion according to frames/comm_sz

  	
 
  unsigned char* pic = new unsigned char[array_size * width * width];
 	 

  //Process 0 will print an error if the number of frames is not a multiple of number of processors.  Exit if true. 
  if (frames%comm_sz!=0){
	if (my_rank==0)
 		printf("Error: #frames not multiple of #processors\n"); 	
        MPI_Finalize();
	return 0;
	}
  if (my_rank==0)
	std::cout<<comm_sz<<"\n";

  // Begginning and end of each array for each processor
  my_start = (array_size*my_rank);
  my_end = (array_size*(my_rank+1));

  MPI_Barrier(MPI_COMM_WORLD);     


  // start time
  struct timeval start, end;
  gettimeofday(&start, NULL);


  // compute frames
  double delta = Delta * pow(0.99,my_start);

  //change the beggining and end that each processor will go through according to the size of their array
  for (int frame = my_start; frame < my_end; frame++) {
    delta *= 0.99;
    const double xMin = xMid - delta;
    const double yMin = yMid - delta;
    const double dw = 2.0 * delta / width;
    for (int row = 0; row < width; row++) {
      const double cy = -yMin - row * dw;
      for (int col = 0; col < width; col++) {
        const double cx = -xMin - col * dw;
        double x = cx;
        double y = cy;
        int depth = 256;
        double x2, y2;
        do {
          x2 = x * x;
          y2 = y * y;
          y = 2 * x * y + cy;
          x = x2 - y2 + cx;
          depth--;
        } while ((depth > 0) && ((x2 + y2) < 5.0));
	
	// (number of frames/number of processors)*rank gives the correct offset 
        pic[(frame-(array_size*my_rank))* width * width + row * width + col] = (unsigned char)depth;  			
      }
    }
  }


  //allocation of new array to be used in gather function to store all values given by each process
  unsigned char* gather_array=NULL;  
if (my_rank==0)
  gather_array= new unsigned char[frames*width*width];
  

  // MPI_Gather function call to collect all data from all processors to process 0
  MPI_Gather(pic, (array_size*width*width), MPI_UNSIGNED_CHAR,gather_array,(array_size*width*width),MPI_UNSIGNED_CHAR, 0,MPI_COMM_WORLD);
  // end time
  gettimeofday(&end, NULL);
  double runtime = end.tv_sec + end.tv_usec / 1000000.0 - start.tv_sec - start.tv_usec / 1000000.0;
  
  if (my_rank==0)
  	printf("compute time: %.4f s\n", runtime);

  // verify result by writing frames to BMP files
  if (my_rank==0){
 	 if ((width <= 400) && (frames <= 30)) {
 	   for (int frame = 0; frame < frames; frame++) {
      		char name[32];
     		 sprintf(name, "fractal%d.bmp", frame + 1000);
     		 writeBMP(width, width, &gather_array[frame * width * width], name);
    }
  }
}
  delete [] pic;

  //finalize all MPI function calls
  MPI_Finalize();

  return 0;
}


