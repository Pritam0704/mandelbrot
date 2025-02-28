#include "bitmap.h"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>

// Function prototypes
int iteration_to_color(int i, int max);
int iterations_at_point(double x, double y, int max);
void *compute_part(void *arg);

// Structure to store thread arguments
typedef struct {
    struct bitmap *bm;
    double xmin, xmax, ymin, ymax;
    int max, start_row, end_row;
} ThreadData;

void show_help()
{
    printf("Use: mandel [options]\n");
    printf("Where options are:\n");
    printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
    printf("-x <coord>  X coordinate of image center point. (default=0)\n");
    printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
    printf("-s <scale>  Scale of the image in Mandelbrot coordinates. (default=4)\n");
    printf("-W <pixels> Width of the image in pixels. (default=500)\n");
    printf("-H <pixels> Height of the image in pixels. (default=500)\n");
    printf("-o <file>   Set output file. (default=mandel.bmp)\n");
    printf("-n <threads> Number of threads to use. (default=1)\n");
    printf("-h          Show this help text.\n\n");
}

int main(int argc, char *argv[])
{
    char c;

    // Default settings
    const char *outfile = "mandel.bmp";
    double xcenter = 0, ycenter = 0, scale = 4;
    int image_width = 500, image_height = 500, max = 1000, num_threads = 1;

    // Parse command-line arguments
    while ((c = getopt(argc, argv, "x:y:s:W:H:m:n:o:h")) != -1) {
        switch (c) {
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
		case 'n': 
			 num_threads = atoi(optarg); 
			 break;	 
        case 'o': 
			outfile = optarg; 
			break;
        
        case 'h':
			 show_help(); 
			 exit(1);
			 break;
        }
    }

    // Display settings
    printf("mandel: x=%lf y=%lf scale=%lf max=%d threads=%d outfile=%s\n",
           xcenter, ycenter, scale, max, num_threads, outfile);

    // Create the image bitmap
    struct bitmap *bm = bitmap_create(image_width, image_height);
    bitmap_reset(bm, MAKE_RGBA(0, 0, 255, 0)); // Default blue background

    // Multi-threading setup
    pthread_t threads[num_threads];
    ThreadData thread_data[num_threads];

    int rows_per_thread = image_height / num_threads;
    
    for (int i = 0; i < num_threads; i++) {
        thread_data[i] = (ThreadData){
            .bm = bm,
            .xmin = xcenter - scale,
            .xmax = xcenter + scale,
            .ymin = ycenter - scale,
            .ymax = ycenter + scale,
            .max = max,
            .start_row = i * rows_per_thread,
            .end_row = (i == num_threads - 1) ? image_height : (i + 1) * rows_per_thread
        };

        pthread_create(&threads[i], NULL, compute_part, &thread_data[i]);
    }

    // Wait for all threads to finish
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    // Save the final image
    if (!bitmap_save(bm, outfile)) {
        fprintf(stderr, "mandel: couldn't write to %s: %s\n", outfile, strerror(errno));
        return 1;
    }

    return 0;
}

/*
 * Compute a part of the Mandelbrot set.
 * Each thread works on a specific set of rows.
 */
void *compute_part(void *arg)
{
    ThreadData *data = (ThreadData *)arg;
    struct bitmap *bm = data->bm;
    double xmin = data->xmin, xmax = data->xmax, ymin = data->ymin, ymax = data->ymax;
    int max = data->max, start = data->start_row, end = data->end_row;

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
 * Compute the number of iterations at a given point in the Mandelbrot set.
 */
int iterations_at_point(double x, double y, int max)
{
    double x0 = x, y0 = y;
    int iter = 0;

    while ((x * x + y * y <= 4) && iter < max) {
        double xt = x * x - y * y + x0;
        double yt = 2 * x * y + y0;
        x = xt;
        y = yt;
        iter++;
    }

    return iteration_to_color(iter, max);
}

/*
 * Convert the number of iterations to a grayscale color.
 */
int iteration_to_color(int i, int max)
{
    int gray = 255 * i / max;
    return MAKE_RGBA(gray, gray, gray, 0);
}
