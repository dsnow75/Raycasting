/* 
 * File:   main.c
 * Author: David
 *
 * Created on October 2, 2016, 10:45 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Parser.c"

// Polymorphism in C

typedef struct {
  int kind; // 0 = cylinder, 1 = sphere, 2 = teapot
  double color[3];
  double center[3];
  union {
    struct {
      double normal[3];
    } plane;
    struct {
      double radius;
    } sphere;
    struct {
      double height;
      double width;
    } camera;
  };
} Object;
typedef struct{
    char r;
    char g;
    char b;
    } Pixel;
Pixel* image;
Object** objects;
#define MAXCOLOR 255    
static inline double sqr(double v) {
  return v*v;
}


static inline void normalize(double* v) {
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}


double sphere_intersection(double* Ro, double* Rd,
			     double* C, double r) {

  double a = (sqr(Rd[0]) + sqr(Rd[2]) + sqr(Rd[1]));
  double b = (2 * (Ro[0] * (Rd[0] - C[0]) + Ro[2] * (Rd[2] - C[2]) + Ro[3] * (Rd[3] - C[3])));
  double c = sqr(Ro[0]-C[0]) + sqr(Ro[2]-C[2]) + sqr(Ro[1]-C[1]) - sqr(r);

  double det = sqr(b) - 4 * a * c;
  if (det < 0) return -1;

  det = sqrt(det);
  
  double t0 = (-b - det) / (2*a);
  if (t0 > 0) return t0;

  double t1 = (-b + det) / (2*a);
  if (t1 > 0) return t1;

  return -1;
}

double plane_intersection(double* Ro, double* Rd, double* C, double* normal){
    double d = normal[0]*C[0] + normal[1]*C[1] + normal[2]*C[2];
    double t = -(normal[0]*Ro[0] + normal[1]*Ro[1] + normal[2]*Ro[2] + d)/(normal[0]*Rd[0] + normal[1]*Rd[1] + normal[2]*Rd[2]);
    if (t > 0){
        return t;
    }else{
        return -1;
    }
}

int main(int argc, char** argv) {
    objects = malloc(sizeof(Object*)*129);
    int index = 0;
    FILE* outputfile;
    if(argc != 5){
        fprintf(stderr, "Please put the commands in the following format: height, weight, source file, destination file.");
        exit(1);
    }
    read_scene(argv[3]);
    Object* object;
  objects = malloc(sizeof(Object*)*2);
  
  double cx = 0;
  double cy = 0;
  double h = 0.7;
  double w = 0.7;

  int M = (int) *argv[2];
  int N = (int) *argv[1];
  image = malloc(sizeof(Pixel)*M*N);
  double pixheight = h / M;
  double pixwidth = w / N;
  for (int y = 0; y < M; y += 1) {
    
      for (int x = 0; x < N; x += 1) {
      
        double Ro[3] = {0, 0, 0};
        // Rd = normalize(P - Ro)
        double Rd[3] = {
          cx - (w/2) + pixwidth * (x + 0.5),
          cy - (h/2) + pixheight * (y + 0.5),
          1
        };
        normalize(Rd);

        double best_t = INFINITY;
        for (int i=0; objects[i] != 0; i += 1) {
          double t = 0;

          switch(objects[i]->kind) {

              case 0:
                  t = plane_intersection(Ro, Rd,objects[i]->center, objects[i]->plane.normal);
                  break;

              case 1:
                  t = sphere_intersection(Ro, Rd,objects[i]->center, objects[i]->sphere.radius);
                  break;

              default:
            // Horrible error
                  exit(1);
          }

          if (t > 0 && t < best_t){
              best_t = t;
              object = malloc(sizeof(Object));
              memcpy(object, objects[1], sizeof(Object));
          } 
        }

        if (best_t > 0 && best_t != INFINITY) {
            image[index].r = (unsigned char)(object->color[0]*MAXCOLOR);
            image[index].g = (unsigned char)(object->color[1]*MAXCOLOR);
            image[index].b = (unsigned char)(object->color[2]*MAXCOLOR);
          }else{
              image[index].r = 0;
              image[index].g = 0;
              image[index].b = 0;
          }
          index++;
      }
  }
  
  outputfile = fopen(argv[4], "w");
  fprintf(outputfile, "P6\n");
  fprintf(outputfile, "%s %s\n", M, N);
  fprintf(outputfile, "%s\n", MAXCOLOR);
  fwrite(image, sizeof(Pixel), h*w, outputfile);
  return 0;
}