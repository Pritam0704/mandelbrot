#include "bitmap.h"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>


int iteration_to_color(int i, int max);
int iterations_at_point(double x, double y, int max);
void *compute_image(void *arg);// This will allow the multiple thread to to execute the funtion

// Structure to store thread arguments 
typedef struct {
    struct bitmap *bm;
    double xmin;
	double xmax;
	double ymin;
	double ymax;
    int max;
	int start_row;
	int end_row;
} 
ThreadData;

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-n <threads> Number of threads to use. (default=1)\n"); // As we gonna add the thread 
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

int main(int argc, char *argv[])
{
    char c;

// These are the default configuration values used
	// if no command line arguments are given.

    const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int    max = 1000;
	int    num_threads = 1;// we initialize the number of threads as 1

    
	// For each command line argument given,
	// override the appropriate configuration value.

    while ((c = getopt(argc, argv, "x:y:s:W:H:m:o:n:h")) != -1)// we added n as the no of the thread
	 {
        switch (c) {
			switch(c) {
				case 'x':
					xcenter = atof(optarg);
					break;
				case 'y':
					ycenter = atof(optarg);
					break;
				case 's':
					scale = atof(optarg);
					break;
				case 'W':
					image_width = atoi(optarg);
					break;
				case 'H':
					image_height = atoi(optarg);
					break;
				case 'm':
					max = atoi(optarg);
					break;
				case 'o':
					outfile = optarg;
					break;
				case 'n': // This will assign the value to num_thread as a integer
					num_threads = atoi(optarg);
					break;
				case 'h':
					show_help();
					exit(1);
					break;
        }
    }
	// Display the configuration of the image.
	// This will print with threads
	printf("mandel: x=%lf y=%lf scale=%lf max=%d threads=%d outfile=%s\n",xcenter,ycenter,scale,max,num_threads,outfile);
    
  // Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

    // Making the Multiple Threads. 

    pthread_t threads[num_threads];
    ThreadData thread_data[num_threads];

    int rows_in_thread = image_height / num_threads;
    
    for (int i = 0; i < num_threads; i++) {
        thread_data[i] = (ThreadData){
            .bm = bm,
            .xmin = xcenter - scale,
            .xmax = xcenter + scale,
            .ymin = ycenter - scale,
            .ymax = ycenter + scale,
            .max = max,
            .start_row = i * rows_in_thread,
            .end_row = (i == num_threads - 1) ? image_height : (i + 1) * rows_in_thread
        };

        pthread_create(&threads[i], NULL, compute_image, &thread_data[i]);
    }

    // waiting to end the threads
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

   // Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}

	return 0;
}

}
/*
 * Mandelbrot set.
 * Each thread works on a specific set of rows.
 * this will run each thread solely and fast.
 */
void *compute_image(void *arg)
{
    ThreadData *data = (ThreadData *)arg;
    struct bitmap *bm = data->bm;
    double xmin = data->xmin;
	double xmax = data->xmax;
	double ymin = data->ymin;
	double ymax = data->ymax;

    int max = data->max;
	int start = data->start_row;
	int end = data->end_row;
    int width = bitmap_width(bm);
    int height = bitmap_height(bm);

    for (int j = start; j < end; j++) {
        for (int i = 0; i < width; i++) {
            double x = xmin + i * (xmax - xmin) / width;
            double y = ymin + j * (ymax - ymin) / height;
            int iters = iterations_at_point(x, y, max);
            bitmap_set(bm, i, j, iters);
        }
    }

    return NULL;
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}

