/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <harrison pettus>
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
// for ease of cout
#include <iostream>
using namespace std;
// added includes
#include <glm/glm.hpp>
#include <vector>


#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;
bool debug = false;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
    glm::vec3 pos;
};

struct Triangle
{
    Vertex v[3];
};

struct Sphere
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
};

struct Light
{
    double position[3];
    double color[3];
};

struct Ray
{
    glm::vec3 position;
    glm::vec3 direction;
    int pixel_W;
    int pixel_H;
    
};

enum objectHit {SPHERE, TRIANGLE};

struct Hit
{
    objectHit type;
    glm::vec3 position;
    float t;
    glm::vec3 normal;
    
    //kd and ks and shin
    glm::vec3 kd;
    glm::vec3 ks;
    float shin;
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];
// added arrays
std::vector<Ray> rays;


int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

glm::vec3 arrayToVec3(double a[3]){
    return glm::vec3(a[0],a[1],a[2]);
}

//generate the rays for each pixel
std::vector<Ray> getRays() {
    //temp vector to store rays created
    std::vector<Ray> temp_rays;
    
    //ray to fill out for each pixel
    //start position is the same for all rays
    Ray r;
    r.position = glm::vec3(0,0,0);
    
    //aspect ratio
    float a = ((float)WIDTH)/((float)HEIGHT);
    if(debug) cout << "w: " << WIDTH << " h: " << HEIGHT << "\na: " << a << endl;
    
    //convert fov to radians
    float fov_radians = ((float)fov) * (M_PI/ 180);
    
    //initial point
    float x = a * tan(fov_radians/2);
    float y = tan(fov_radians/2);
    float z = -1.0f;
    if(debug) cout << "x: " << x << "\ty: " << y << endl;
    
    for(size_t i = 0; i<WIDTH; i++) {
        for(size_t j = 0; j<HEIGHT; j++) {
            glm::vec3 temp;
            temp.x = -x + (i / (float)WIDTH) * 2 * x;
            temp.y = -y + (j / (float)HEIGHT) * 2 * y;
            temp.z = -1; //this stays const
            
            r.direction = glm::normalize(temp);
            r.pixel_W = i;
            r.pixel_H = j;
            // add ray at location i,j to ray vector
            temp_rays.push_back(r);
        }
    }
    
    return temp_rays;
}

float hit_sphere(Ray& r, Sphere& s, Hit& hit) {
    bool negate = false;
    
    glm::vec3 center = arrayToVec3(s.position);
    double radius = s.radius;
    glm::vec3 oc = r.position - center;
    float a = dot(r.direction, r.direction);
    float b = 2.0 * glm::dot(oc, r.direction);
    float c = glm::dot(oc, oc) - (radius*radius);
    float discriminant = b*b - 4*c;
    float hit_loc;
    if(discriminant > 0) {
        float t0 = -b + sqrt(discriminant);
        t0 /=2.0f;
        float t1 = -b - sqrt(discriminant);
        t1 /=2.0f;
        hit_loc = min(t0,t1);
        if (t0 < 0 && t1 <0) return -1;
        if (t0 < 0 && t1 >=0) {
            hit_loc = t1;
            //cout << "here";
            negate = true;
        }
        if (t0 >= 0 && t1 <0) {
            hit_loc = t0;
            //cout << "here";
            negate = true;
        }
        //cout << t0 << " " << t1 << "\tmin_t = " << hit_loc << endl;
    }
    else
        return -1;
    
    
    //fill out hit info
    hit.type = SPHERE;
    hit.position = r.position + r.direction * hit_loc;
    hit.t = hit_loc;
    hit.normal = glm::normalize(hit.position - center);
    if(negate) hit.normal *= -1;
    
    hit.kd = arrayToVec3(s.color_diffuse);
    //cout << hit.kd.x << " " << hit.kd.y << " " << hit.kd.z << endl;
    hit.ks = arrayToVec3(s.color_specular);
    //cout << hit.ks.x << " " << hit.ks.y << " " << hit.ks.z << endl;
    hit.shin = s.shininess;
    
    
    return hit_loc;
}

float hit_triangle(Ray& r, Triangle &t, Hit &hit) {
    glm::vec3 v0 = arrayToVec3(t.v[0].position);
    glm::vec3 v1 = arrayToVec3(t.v[1].position);
    glm::vec3 v2 = arrayToVec3(t.v[2].position);
    
    // plane normal containing triangle
    glm::vec3 v0v1 = v1 - v0;
    glm::vec3 v0v2 = v2 - v0;
    glm::vec3 normal = glm::cross(v0v1, v0v2);
    float area2 = glm::length(normal);
    normal = glm::normalize(normal);
    
    float nDotD = glm::dot(normal, r.direction);
    
    //they don't intersect
    if(fabs(nDotD) < 1e-6){
        return -1;
    }
    float d = -glm::dot(normal, v0);
    float tt = -(glm::dot(normal, r.position)+d) / glm::dot(normal, r.direction);
    //hit.position = r.position + r.direction * hit_loc;
    if (tt <= 0) return -1; // the triangle is behind
    
    // point of intersect plane
    glm::vec3 intersect_point = r.position + tt * r.direction;

    
    
    glm::vec3 edge0 = v1 - v0;
    glm::vec3 edge1 = v2 - v0;
    glm::vec3 edge2 = intersect_point - v0;
    float d00 = glm::dot(edge0, edge0);
    float d01 = glm::dot(edge0, edge1);
    float d11 = glm::dot(edge1, edge1);
    float d20 = glm::dot(edge2, edge0);
    float d21 = glm::dot(edge2, edge1);
    float denom = (d00 * d11) - (d01 * d01);
    if (fabs(denom) < 1e-6) {
        //intersect outside traingle
        return -1;
    }
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1 - v - w;
    if (u < 0 || v < 0 || w < 0 || u > 1 || v > 1 || w > 1) {
        //intersect outside traingle
        return -1;
    }
    
    
    //fill out Hit object
    hit.type = TRIANGLE;
    hit.position = intersect_point;
    hit.t = tt;
    hit.normal = glm::normalize(u * arrayToVec3(t.v[0].normal) + v * arrayToVec3(t.v[1].normal) + w * arrayToVec3(t.v[2].normal));
    hit.kd = u * arrayToVec3(t.v[0].color_diffuse) + v * arrayToVec3(t.v[1].color_diffuse) + w * arrayToVec3(t.v[2].color_diffuse);
    if (hit.kd != glm::vec3(0,0,0)) {
        hit.kd = glm::clamp(hit.kd, 0.0f, 1.0f);
    }
    
    //cout << hit.kd.x << " " << hit.kd.y << " " << hit.kd.z << endl;
    hit.ks = (u * arrayToVec3(t.v[0].color_specular) + v * arrayToVec3(t.v[1].color_specular) + w * arrayToVec3(t.v[2].color_specular));
    
    if (hit.ks != glm::vec3(0,0,0)) {
        hit.ks = glm::clamp(hit.ks, 0.0f, 1.0f);
    }
    hit.shin = u * t.v[0].shininess + v * t.v[1].shininess + w * t.v[2].shininess;
    
    
    return tt;
    
}


Hit closestHit(Ray ray) {
    std::vector<Hit> allHits;
    Hit closestHit;
    
    for(size_t i=0; i<num_triangles; i++) {
        Hit hit2;
        // get t location of hit
        float t = hit_triangle(ray, triangles[i], hit2);
        // there was a hit need to check distance to confirm
        if(t != -1){
            // get vector to hit
            glm::vec3 pq = ray.position - hit2.position;
            //avoiding spotting
            if(glm::distance(ray.position, hit2.position) > .0002) {
                allHits.push_back(hit2);
            }
        }
    }
    
    for(size_t i=0; i<num_spheres; i++) {
        Hit hit;
        // get t location of hit
        float t = hit_sphere(ray, spheres[i], hit);
        
        //if it did hit the sphere
        if(t != -1){
            // get vector to hit
            glm::vec3 pq = ray.position - hit.position;
            //avoiding spotting
            if(glm::distance(ray.position, hit.position) > .0008) {
                allHits.push_back(hit);
            }
        }
    }
    
    
    
    if (allHits.size() > 0) {
        closestHit.t = 10000;
        for(size_t i=0; i < allHits.size(); i++) {
            if (allHits[i].t < closestHit.t) {
                
                closestHit = allHits[i];
            }
        }
        return closestHit;
    }
    
    closestHit.t = -100;
    return closestHit;
    
    
}


glm::vec3 getColor(Hit hit) {
    
    //final color (summation from all lights)
    glm::vec3 fin = glm::vec3(0,0,0);
    
    for(size_t x=0; x<num_lights; x++) {
        glm::vec3 lightSource = arrayToVec3(lights[x].position);
        glm::vec3 lightColor = arrayToVec3(lights[x].color);
        glm::vec3 n = hit.normal;
        glm::vec3 l = glm::normalize(lightSource - hit.position);
        glm::vec3 v = glm::normalize(-hit.position);
        glm::vec3 r = 2 * (glm::dot(l, n)) * n - l;
        //glm::vec3 r = glm::reflect(l,n);
        r = glm::normalize(r);
        
        //clamp values for division
        float lDotN = glm::dot(l,n);
        if (lDotN < 0) lDotN = 0;
        float rDotV = glm::dot(r,v);
        if (rDotV < 0) rDotV = 0;
        
        //ray from hit position to light source
        Ray shadowRay;
        shadowRay.position = hit.position;
        shadowRay.direction = l;
        
        //set upper bound for t
        float max_t = (arrayToVec3(lights[x].position).y - shadowRay.position.y) / shadowRay.direction.y;
        
        // shadows
        Hit p = closestHit(shadowRay);
        
        if (p.t != -100 && (p.t < max_t)) {
            //cout << p.t << endl;
            continue;
        }
        else {
            glm::vec3 col = lightColor * (hit.kd * (lDotN) + hit.ks * pow(rDotV, hit.shin));
            //col = glm::normalize(col);
            fin += (col);
        }
        
    }
    
    return fin;
}

//MODIFY THIS FUNCTION
void draw_scene()
{
    //set generated rays to the rays vector
    rays = getRays();
    if(debug) cout << rays.size();
    
    
    int rayCount =0;
    //a simple test output
    for(unsigned int x=0; x<WIDTH; x++)
    {
        glPointSize(2.0);
        glBegin(GL_POINTS);
        for(unsigned int y=0; y<HEIGHT; y++)
        {
            Hit closestIntPoint = closestHit(rays[rayCount]);
            
            //did intersect with object
            if (closestIntPoint.t != -100) {
                rayCount++;
                glm::vec3 finalColor = getColor(closestIntPoint);

                //add the ambient light once per pixel
                finalColor += arrayToVec3(ambient_light);
                finalColor = glm::clamp(finalColor, 0.0f ,1.0f);
                plot_pixel(x, y, finalColor.x*255.0f, finalColor.y*255.0f, finalColor.z*255.0f);
            }
            
            // no collisions, color background color
            else {
                rayCount++;
                plot_pixel(x, y, 255, 255, 255);
            }
            
        }
        glEnd();
        glFlush();
    }
    printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
      glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
      glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
      buffer[y][x][0] = r;
      buffer[y][x][1] = g;
      buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    plot_pixel_display(x,y,r,g,b);
    if(mode == MODE_JPEG)
        plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
    printf("Saving JPEG file: %s\n", filename);

    ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
    if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
        printf("Error in Saving\n");
    else
        printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
    if(strcasecmp(expected,found))
    {
        printf("Expected '%s ' found '%s '\n", expected, found);
        printf("Parse error, abnormal abortion\n");
        exit(0);
    }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check(check,str);
    fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
    printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check("rad:",str);
    fscanf(file,"%lf",r);
    printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
    char s[100];
    fscanf(file,"%s",s);
    parse_check("shi:",s);
    fscanf(file,"%lf",shi);
    printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
    FILE * file = fopen(argv,"r");
    int number_of_objects;
    char type[50];
    Triangle t;
    Sphere s;
    Light l;
    fscanf(file,"%i", &number_of_objects);

    printf("number of objects: %i\n",number_of_objects);

    parse_doubles(file,"amb:",ambient_light);

    for(int i=0; i<number_of_objects; i++)
    {
        fscanf(file,"%s\n",type);
        printf("%s\n",type);
        if(strcasecmp(type,"triangle")==0)
        {
            printf("found triangle\n");
            for(int j=0;j < 3;j++)
            {
                parse_doubles(file,"pos:",t.v[j].position);
                parse_doubles(file,"nor:",t.v[j].normal);
                parse_doubles(file,"dif:",t.v[j].color_diffuse);
                parse_doubles(file,"spe:",t.v[j].color_specular);
                parse_shi(file,&t.v[j].shininess);
            }

            if(num_triangles == MAX_TRIANGLES)
            {
                printf("too many triangles, you should increase MAX_TRIANGLES!\n");
                exit(0);
            }
            triangles[num_triangles++] = t;
        }
        else if(strcasecmp(type,"sphere")==0)
        {
            printf("found sphere\n");

            parse_doubles(file,"pos:",s.position);
            parse_rad(file,&s.radius);
            parse_doubles(file,"dif:",s.color_diffuse);
            parse_doubles(file,"spe:",s.color_specular);
            parse_shi(file,&s.shininess);

            if(num_spheres == MAX_SPHERES)
            {
                printf("too many spheres, you should increase MAX_SPHERES!\n");
                exit(0);
            }
            spheres[num_spheres++] = s;
        }
        else if(strcasecmp(type,"light")==0)
        {
            printf("found light\n");
            parse_doubles(file,"pos:",l.position);
            parse_doubles(file,"col:",l.color);

            if(num_lights == MAX_LIGHTS)
            {
                printf("too many lights, you should increase MAX_LIGHTS!\n");
                exit(0);
            }
            lights[num_lights++] = l;
        }
        else
        {
            printf("unknown type in scene description:\n%s\n",type);
            exit(0);
        }
    }
    return 0;
}

void display()
{
    
}

void init()
{
    glMatrixMode(GL_PROJECTION);
    glOrtho(0,WIDTH,0,HEIGHT,1,-1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
    //hack to make it only draw once
    static int once=0;
    if(!once)
    {
        draw_scene();
        if(mode == MODE_JPEG)
            save_jpg();
    }
    once=1;
}

int main(int argc, char ** argv)
{
    if ((argc < 2) || (argc > 3))
    {
        printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
        exit(0);
    }
    if(argc == 3)
    {
        mode = MODE_JPEG;
        filename = argv[2];
    }
    else if(argc == 2)
        mode = MODE_DISPLAY;

    glutInit(&argc,argv);
    loadScene(argv[1]);

    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(WIDTH,HEIGHT);
    int window = glutCreateWindow("Ray Tracer");
    #ifdef __APPLE__
        // This is needed on recent Mac OS X versions to correctly display the window.
        glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
    #endif
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    init();
    glutMainLoop();
}

