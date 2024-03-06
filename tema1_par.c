// Author: APD team, except where source was noted

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "helpers.h"

#define CONTOUR_CONFIG_COUNT 16
#define FILENAME_MAX_SIZE 50
#define STEP 8
#define SIGMA 200
#define RESCALE_X 2048
#define RESCALE_Y 2048

#define CLAMP(v, min, max) \
  if (v < min) {           \
    v = min;               \
  } else if (v > max) {    \
    v = max;               \
  }

typedef struct
{
  int id;
  int NR_THREADS;
  ppm_image *image;
  ppm_image *scaled_image;
  ppm_image **contour_map;
  unsigned char **grid;
  pthread_barrier_t *bariera;
  uint8_t sample[3];
  int step_x;
  int step_y;

} threads_aux_info;

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
  ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
  if (!map) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }

  for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
    char filename[FILENAME_MAX_SIZE];
    sprintf(filename, "./contours/%d.ppm", i);
    map[i] = read_ppm(filename);
  }

  return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
  for (int i = 0; i < contour->x; i++) {
    for (int j = 0; j < contour->y; j++) {
      int contour_pixel_index = contour->x * i + j;
      int image_pixel_index = (x + i) * image->y + y + j;

      image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
      image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
      image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
    }
  }
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
unsigned char **sample_grid(ppm_image *image, int step_x, int step_y, unsigned char sigma) {
  int p = image->x / step_x;
  int q = image->y / step_y;

  unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char *));
  if (!grid) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }

  for (int i = 0; i <= p; i++) {
    grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
    if (!grid[i]) {
      fprintf(stderr, "Unable to allocate memory\n");
      exit(1);
    }
  }

  for (int i = 0; i < p; i++) {
    for (int j = 0; j < q; j++) {
      ppm_pixel curr_pixel = image->data[i * step_x * image->y + j * step_y];

      unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

      if (curr_color > sigma) {
        grid[i][j] = 0;
      } else {
        grid[i][j] = 1;
      }
    }
  }
  grid[p][q] = 0;

  // last sample points have no neighbors below / to the right, so we use pixels on the
  // last row / column of the input image for them
  for (int i = 0; i < p; i++) {
    ppm_pixel curr_pixel = image->data[i * step_x * image->y + image->x - 1];

    unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

    if (curr_color > sigma) {
      grid[i][q] = 0;
    } else {
      grid[i][q] = 1;
    }
  }
  for (int j = 0; j < q; j++) {
    ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * step_y];

    unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

    if (curr_color > sigma) {
      grid[p][j] = 0;
    } else {
      grid[p][j] = 1;
    }
  }

  return grid;
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(ppm_image *image, unsigned char **grid, ppm_image **contour_map, int step_x, int step_y) {
  int p = image->x / step_x;
  int q = image->y / step_y;

  for (int i = 0; i < p; i++) {
    for (int j = 0; j < q; j++) {
      unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
      update_image(image, contour_map[k], i * step_x, j * step_y);
    }
  }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
  for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
    free(contour_map[i]->data);
    free(contour_map[i]);
  }
  free(contour_map);

  for (int i = 0; i <= image->x / step_x; i++) {
    free(grid[i]);
  }
  free(grid);

  free(image->data);
  free(image);
}



ppm_image *rescale_image(ppm_image *image) {
  uint8_t sample[3];

  // we only rescale downwards
  if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
    return image;
  }

  // alloc memory for image
  ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));
  if (!new_image) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }
  new_image->x = RESCALE_X;
  new_image->y = RESCALE_Y;

  new_image->data = (ppm_pixel *)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
  if (!new_image) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }

  // use bicubic interpolation for scaling
  for (int i = 0; i < new_image->x; i++) {
    for (int j = 0; j < new_image->y; j++) {
      float u = (float)i / (float)(new_image->x - 1);
      float v = (float)j / (float)(new_image->y - 1);
      sample_bicubic(image, u, v, sample);

      new_image->data[i * new_image->y + j].red = sample[0];
      new_image->data[i * new_image->y + j].green = sample[1];
      new_image->data[i * new_image->y + j].blue = sample[2];
    }
  }

  free(image->data);
  free(image);

  return new_image;
}

// Alocates memory for the ppm_image structure 
ppm_image *create_rescaled_image() {
  ppm_image *image = (ppm_image *)malloc(sizeof(ppm_image));
  if (!image) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }
  image->x = RESCALE_X;
  image->y = RESCALE_Y;
  image->data = (ppm_pixel *)malloc(image->x * image->y * sizeof(ppm_pixel));
  image->data = (ppm_pixel *)malloc(image->x * image->y * sizeof(ppm_pixel));
  if (!image) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }

    return image;
}

// Alocates memory for the grid 
unsigned char **create_grid(int p, int q){
  unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char *));
  if (!grid) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }

  for (int i = 0; i <= p; i++) {
    grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
    if (!grid[i]) {
      fprintf(stderr, "Unable to allocate memory\n");
      exit(1);
    }
  }
    return grid;

}

void *thread_function(void *args) {
  threads_aux_info info = *(threads_aux_info *)args;
  int start_rescale = info.id * ceil((double)info.scaled_image->x / info.NR_THREADS);
  int end_rescale = fmin((info.id + 1) * ceil((double)info.scaled_image->x / info.NR_THREADS), info.scaled_image->x);

  // rescale image
  if (info.image->x <= RESCALE_X && info.image->y <= RESCALE_Y) {
    info.scaled_image = info.image;
  } else {
    for (int i = start_rescale; i < end_rescale; i++) {
      for (int j = 0; j < info.scaled_image->y; j++) {
        float u = (float)i / (float)(info.scaled_image->x - 1);
        float v = (float)j / (float)(info.scaled_image->y - 1);
        sample_bicubic(info.image, u, v, info.sample);

        info.scaled_image->data[i * info.scaled_image->y + j].red = info.sample[0];
        info.scaled_image->data[i * info.scaled_image->y + j].green = info.sample[1];
        info.scaled_image->data[i * info.scaled_image->y + j].blue = info.sample[2];
      }
    }
    // wait for the whole image to be rescaled
    pthread_barrier_wait(info.bariera);
  }

  // sample the grid
  int p = info.scaled_image->x / info.step_x;
  int q = info.scaled_image->y / info.step_y;

  int start_p = info.id * ceil((double)p / info.NR_THREADS);
  int end_p = fmin((info.id + 1) * ceil((double)p / info.NR_THREADS), p);

  int start_q = info.id * ceil((double)p / info.NR_THREADS);
  int end_q = fmin((info.id + 1) * ceil((double)q / info.NR_THREADS), q);

  for (int i = start_p; i < end_p; i++) {
    for (int j = 0; j < q; j++) {
      ppm_pixel curr_pixel = info.scaled_image->data[i * info.step_x * info.scaled_image->y + j * info.step_y];
      unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

      if (curr_color > SIGMA) {
        info.grid[i][j] = 0;
      } else {
        info.grid[i][j] = 1;
      }
    }
  }
  info.grid[p][q] = 0;

  // last sample points have no neighbors below / to the right, so we use pixels on the
  // last row / column of the input image for them
  for (int i = start_p; i < end_p; i++) {
    ppm_pixel curr_pixel = info.scaled_image->data[i * info.step_x * info.scaled_image->y + info.scaled_image->x - 1];
    unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

    if (curr_color > SIGMA) {
      info.grid[i][q] = 0;
    } else {
      info.grid[i][q] = 1;
    }
  }
  for (int j = start_q; j < end_q; j++) {
    ppm_pixel curr_pixel = info.scaled_image->data[(info.scaled_image->x - 1) * info.scaled_image->y + j * info.step_y];
    unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

    if (curr_color > SIGMA) {
      info.grid[p][j] = 0;
    } else {
      info.grid[p][j] = 1;
    }
  }

  // wait for the whole grid to be created
  pthread_barrier_wait(info.bariera);

  // march the squares

  for (int i = start_p; i < end_p; i++) {
    for (int j = 0; j < q; j++) {
      unsigned char k = 8 * info.grid[i][j] + 4 * info.grid[i][j + 1] + 2 * info.grid[i + 1][j + 1] + 1 * info.grid[i + 1][j];
      update_image(info.scaled_image, info.contour_map[k], i * info.step_x, j * info.step_y);
    }
  }

  pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
  if (argc < 4) {
    fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
    return 1;
  }

  int NR_OF_THREADS = atoi(argv[3]);
  pthread_t threads[NR_OF_THREADS];
  threads_aux_info threads_info[NR_OF_THREADS];
  int id_thread;
  int err;
  void *status = NULL;
  pthread_barrier_t bariera_threaduri;
  pthread_barrier_init(&bariera_threaduri, NULL, NR_OF_THREADS);

  ppm_image **contour_map = init_contour_map();
  ppm_image *image = read_ppm(argv[1]);
  ppm_image *scaled_image = NULL;
  int step_x = STEP;
  int step_y = STEP;

  if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
    scaled_image = image;
  } else {
    scaled_image = create_rescaled_image();
  }

  int p = scaled_image->x / step_x;
  int q = scaled_image->y / step_y;
  unsigned char **grid = create_grid(p,q);

  // create threads
  for (id_thread = 0; id_thread < NR_OF_THREADS; id_thread++) {
    threads_info[id_thread].id = id_thread;
    threads_info[id_thread].bariera = &bariera_threaduri;
    threads_info[id_thread].scaled_image = scaled_image;
    threads_info[id_thread].contour_map = contour_map;
    threads_info[id_thread].image = image;
    threads_info[id_thread].grid = grid;
    threads_info[id_thread].NR_THREADS = NR_OF_THREADS;
    threads_info[id_thread].step_x = step_x;
    threads_info[id_thread].step_y = step_y;

    err = pthread_create(&threads[id_thread], NULL, thread_function, &threads_info[id_thread]);

    if (err) {
      printf("Error creating the thread %d\n", id_thread);
      exit(-1);
    }
  }

  // join threads
  for (id_thread = 0; id_thread < NR_OF_THREADS; id_thread++) {
    err = pthread_join(threads[id_thread], &status);

    if (err) {
      printf("Error waiting for thread to finish  %d\n", id_thread);
      exit(-1);
    }
  }

  pthread_barrier_destroy(&bariera_threaduri);

  // write output
  write_ppm(scaled_image, argv[2]);

  
  if (scaled_image != image) {
    free(image->data);
    free(image);
  }
  
  free_resources(scaled_image, contour_map, grid, step_x);

  return 0;
}
