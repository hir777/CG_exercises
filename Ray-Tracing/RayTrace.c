#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "Scene.h"
#include "RayTrace.h"
#include "Geometry.h"

#define MIN(a,b) (((a)<(b))?(a):(b))


static int rayTrace(Vector3, Vector3,
		    Scene, Color*);


// Initialize the image with the background color specified as input.
// width and height corresponds to the dimension of the image.
static void
initImageWithBackground(Color background_color, Color** image,
			int width, int height)
{
  int i;
  int j;
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      image[i][j]._red = background_color._red;
      image[i][j]._green = background_color._green;
      image[i][j]._blue = background_color._blue;
    }
  }
}


// Clamp c's entries between m and M. 
static void clamp(Color* c, float m, float M) {
  c->_red = fminf(fmaxf(c->_red, m), M);
  c->_green = fminf(fmaxf(c->_green, m), M);
  c->_blue = fminf(fmaxf(c->_blue, m), M);
}


// Complete
// Given a ray (origin, direction), check if it intersects a given
// sphere.
// Return 1 if there is an intersection, 0 otherwise.
// *t contains the distance to the closest intersection point, if any.
static int hitSphere(Vector3 origin, Vector3 direction, Sphere sphere, float* t) {
  float a, c, b_square, d;
  float dPlus, dMinus, dist;
  Vector3 tmp;

  normalize(direction, &direction);

  // a 
  sub(sphere._center, origin, &tmp);
  computeDotProduct(direction, tmp, &a);

  // c
  sub(origin, sphere._center, &tmp);
  computeNorm(tmp, &c);

  // b ^ 2
  b_square = powf(c, 2.0) - powf(a, 2.0);

  // Check if there is an intersection.
  if( powf(sphere._radius, 2.0) - b_square < 0.0 )
    return 0;                                     // Return 0 only when there is no intersection.

  
  dPlus = sqrtf( powf(sphere._radius, 2.0) - b_square );
  dMinus = -1.0 * dPlus;

  if( a - dPlus > 0 )
    dist = MIN( a - dPlus, a - dMinus);
  else
    dist = a - dMinus;

  mulAV(dist, direction, &tmp);

  add(origin, tmp, &tmp);

  computeNorm(tmp, t);

  return 1;
}


// Check if the ray defined by (scene._camera, direction) is intersecting
// any of the spheres defined in the scene.
// Return 0 if there is no intersection, and 1 otherwise.
//
// If there is an intersection:
// - the position of the intersection with the closest sphere is computed 
// in hit_pos
// - the normal to the surface at the intersection point in hit_normal
// - the diffuse color and specular color of the intersected sphere
// in hit_color and hit_spec
static int
hitScene(Vector3 origin, Vector3 direction, Scene scene,
	 Vector3* hit_pos, Vector3* hit_normal,
	 Color* hit_color, Color* hit_spec)
{
  Vector3 o = origin;
  Vector3 d = direction;

  float t_min = FLT_MAX;
  int hit_idx = -1;
  Sphere hit_sph;

  // For each sphere in the scene
  int i;
  for (i = 0; i < scene._number_spheres; ++i) {
    Sphere curr = scene._spheres[i];
    float t = 0.0f;
    if (hitSphere(o, d, curr, &t)) {
      if (t < t_min) {
	hit_idx = i;
	t_min = t;
	hit_sph = curr;
      }
    }
  }

  if (hit_idx == -1) return 0;

  Vector3 td;
  mulAV(t_min, d, &td);
  add(o, td, hit_pos);
  
  Vector3 n;
  sub(*hit_pos, hit_sph._center, &n);
  mulAV(1.0f / hit_sph._radius, n, hit_normal);

  // Save the color of the intersected sphere in hit_color and hit_spec
  *hit_color = hit_sph._color;
  *hit_spec = hit_sph._color_spec;
  
  return 1;
}


// Save the image in a raw buffer (texture)
// The memory for texture is allocated in this function. It needs to 
// be freed in the caller.
static void saveRaw(Color** image, int width, int height, GLubyte** texture) {
  int count = 0;
  int i;
  int j;
  *texture = (GLubyte*)malloc(sizeof(GLubyte) * 3 * width * height);
	
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      unsigned char red = (unsigned char)(image[i][j]._red * 255.0f);
      unsigned char green = (unsigned char)(image[i][j]._green * 255.0f);
      unsigned char blue = (unsigned char)(image[i][j]._blue * 255.0f);
			
      (*texture)[count] = red;
      count++;
			
      (*texture)[count] = green;
      count++;
			
      (*texture)[count] = blue;
      count++;
    }
  }
}


// Complete
// Given an intersection point (hit_pos),
// the normal to the surface at the intersection point (hit_normal),
// and the color (diffuse and specular) terms at the intersection point,
// compute the color intensity at the point by applying the Phong
// shading model.
// Return the color intensity in *color.
static void shade(Vector3 hit_pos, Vector3 hit_normal,
      Color hit_color, Color hit_spec, Scene scene, Color* color)
{
  Vector3 L, V, R, I, toLight;
  Vector3 pos, normal;
  Color cl, spec;
  float t;

  // Complete
  // ambient component

  color->_red   += hit_color._red * scene._ambient._red;
  color->_green += hit_color._green * scene._ambient._green;
  color->_blue  += hit_color._blue * scene._ambient._blue; 

   // for each light in the scene
  int l;
  for (l = 0; l < scene._number_lights; l++) {
    // Complete
    // Form a shadow ray and check if the hit point is under
    // direct illumination from the light source

    int i, hasObstacle = 0;
    sub(scene._lights[l]._light_pos, hit_pos, &toLight);
    normalize(toLight, &toLight);
    hasObstacle = hitScene(hit_pos, toLight, scene, &pos, &normal, &cl, &spec);

    // Check if there is an obstacle in the direction from hit position to the light source
    if( hasObstacle == 0 ) {

        // Complete
        // diffuse component

        // Get the unit vector L which is from current point to the light source.
        sub(scene._lights[l]._light_pos, hit_pos, &L);
        normalize(L, &L);

        // Calculate Id = kd * I * max(N・L, 0)
        computeDotProduct(hit_normal, L, &t);
        t = fmaxf( t, 0.0);
        color->_red   += hit_color._red * scene._lights[l]._light_color._red * t;
        color->_green += hit_color._green * scene._lights[l]._light_color._green * t;
        color->_blue  += hit_color._blue * scene._lights[l]._light_color._blue * t;


        // Complete
        // specular component    Phong-model

        // Get the unit vector V which is the viewing direction
        sub(scene._camera, hit_pos, &V);
        normalize(V, &V);

        // Get the unit vector I from the light source point to the hit position
        sub(hit_pos, scene._lights[l]._light_pos, &I);
        normalize(I, &I);

        // Get the unit vector directed towards the ideal specular reflection
        computeDotProduct(I, hit_normal, &t);
        t *= -2.0;
        mulVA(hit_normal, t, &R);
        add(I, R, &R);

        computeDotProduct(V, R, &t);
        t = pow( fmaxf(t, 0.0), 10);         // You have to fine-tune the shininess value

        color->_red   += hit_spec._red * scene._lights[l]._light_color._red * t;
        color->_green += hit_spec._green * scene._lights[l]._light_color._green * t;
        color->_blue  += hit_spec._blue * scene._lights[l]._light_color._blue * t;

    }

  }
}


static int rayTrace(Vector3 origin, Vector3 direction_normalized,
		    Scene scene, Color* color)
{
  Vector3 hit_pos;
  Vector3 hit_normal;
  Color hit_color;
  Color hit_spec;
  int hit;
    
  // does the ray intersect an object in the scene?
  hit = 
    hitScene(origin, direction_normalized, scene,
	     &hit_pos, &hit_normal, &hit_color,
	     &hit_spec);
  
  // no hit
  if (!hit)
    return 0;
  color->_red = 0;
  color->_green = 0;
  color->_blue = 0;

  // otherwise, apply the shading model at the intersection point
  shade(hit_pos, hit_normal, hit_color, hit_spec, scene, color);

  return 1;
}


void rayTraceScene(Scene scene, int width, int height, GLubyte** texture) {
  Color** image;
  int i;
  int j;
	
  Vector3 camera_pos;
  float screen_scale;
  
  image = (Color**)malloc(height * sizeof(Color*));
  for (i = 0; i < height; i++) {
    image[i] = (Color*)malloc(width * sizeof(Color));
  }
	
  // Init the image with the background color
  initImageWithBackground(scene._background_color, image, width, height);
  
  // get parameters for the camera position and the screen fov
  camera_pos._x = scene._camera._x;
  camera_pos._y = scene._camera._y;
  camera_pos._z = scene._camera._z;
	
  screen_scale = scene._scale;
	
	
  // go through each pixel
  // and check for intersection between the ray and the scene
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      // Compute (x,y) coordinates for the current pixel 
      // in scene space
      float x = screen_scale * j - 0.5f * screen_scale * width;
      float y = screen_scale * i - 0.5f * screen_scale * height;

      // Form the vector camera to current pixel
      Vector3 direction;
      Vector3 direction_normalized;
      
      direction._x = x - camera_pos._x;
      direction._y = y - camera_pos._y;
      direction._z = -camera_pos._z;
			
      normalize(direction, &direction_normalized);

      Vector3 origin = scene._camera;
      Color color;
      int hit = rayTrace(origin, direction_normalized, scene, &color);
      if (hit) {
	      clamp(&color, 0.f, 1.f);
	      image[i][j] = color;
      }

    }

  }

  // save image to texture buffer
  saveRaw(image, width, height, texture);

  for (i = 0; i < height; i++) {
    free(image[i]);
  }
	
  free(image);
}
