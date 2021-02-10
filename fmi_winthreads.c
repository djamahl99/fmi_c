#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_resize.h"
#include <stdbool.h>
#include <stdio.h>
// #include <curl/curl.h>
// #include <curl/types.h>
// #include <curl/easy.h>
#include <string.h>
#include <windows.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <conio.h>
#include <process.h>
#define DESIRED_COLOURS 50
#define MAX_SIZE 2500

void ThreadFunc(void* args);

HANDLE hColoursMutex;

double* RGBtoXYZ(unsigned char* rgb) {
    double r = rgb[0] / 255.0;
    double g = rgb[1] / 255.0;
    double b = rgb[2] / 255.0;

    if (r > 0.04045) {
        r = pow(((r + 0.055) / 1.055), 2.4);
    } else {
        r = r / 12.92;
    }

    if (g > 0.04045) {
        g = pow(((g + 0.055) / 1.055), 2.4);
    } else {
        g = g / 12.92;
    }

    if (b > 0.04045) {
        b = pow(((b + 0.055) / 1.055), 2.4);
    } else {
        b = b / 12.92;
    }

    r *= 100;
    g *= 100;
    b *= 100;

    double* ret = calloc(3, sizeof(double));
    ret[0] = r * 0.4124 + g * 0.3576 + b * 0.1805;
    ret[1] = r * 0.2126 + g * 0.7152 + b * 0.0722;
    ret[2] = r * 0.0193 + g * 0.1192 + b * 0.9505;
    return ret;
}

double* XYZtoCIELab(double* xyz) {
	double refX = 95.047;
    double refY = 100;
    double refZ = 108.883;
    double x = xyz[0] / refX;
    double y = xyz[1] / refY;
    double z = xyz[2] / refZ;


    if (x > 0.008856) {
        x = pow(x, 0.333333333333);

    } else {
        x = (7.787 * x) + (16 / 116);
    }

    if (y > 0.008856) {
        y = pow(y, 0.333333333333);
    } else {
        y = (7.787 * y) + (16 / 116);
    }

    if (z > 0.008856) {
        z = pow(z, 0.333333333333);
    } else {
        z = (7.787 * z) + (16 / 116);
    }


    double* ret = calloc(3, sizeof(double));
    ret[0] = (116 * y) - 16;
    ret[1] = 500 * (x - y);
    ret[2] = 200 * (y - z);
    return ret;
}

double labEuclideanDistance(unsigned char* c1, unsigned char* c2) {
    double* xyz1 = RGBtoXYZ(c1);
    double* xyz2 = RGBtoXYZ(c2);
    double* lab1 = XYZtoCIELab(xyz1);
    double* lab2 = XYZtoCIELab(xyz2);

    double dis = sqrt(pow(lab2[0] - lab1[0], 2) + pow(lab2[1] - lab1[1], 2) + pow(lab2[2] - lab1[2], 2));

    free(xyz1);
    free(xyz2);
    free(lab1);
    free(lab2);

    return dis;
}

double labDist(unsigned char r1, unsigned char g1, unsigned char b1, unsigned char r2, unsigned char g2, unsigned char b2) {
    unsigned char c1[3];
    c1[0] = r1;
    c1[1] = g1;
    c1[2] = b1;

    unsigned char c2[3];
    c2[0] = r2;
    c2[1] = g2;
    c2[2] = b2;

    return labEuclideanDistance(c1, c2);
}

double labComp(unsigned char* c1, unsigned char* c2) {
    double* xyz1 = RGBtoXYZ(c1);
    double* xyz2 = RGBtoXYZ(c2);
    double* lab1 = XYZtoCIELab(xyz1);
    double* lab2 = XYZtoCIELab(xyz2);

    double dis = (lab2[0] - lab1[0]) + (lab2[1] - lab1[1]) + (lab2[2] - lab1[2]);
    
    free(xyz1);
    free(xyz2);
    free(lab1);
    free(lab2);

    return dis;
}

typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
    int count;
    int trueIndex;
    int colorIndex;
    double metricDiff;
} Colour;

unsigned int getIndex_(Colour* c) {
    return getIndex(c->r, c->g, c->b);
}


double labDistStruct(Colour* c1, Colour* c2) {
    // printf("c1 %d %d %d , c2 %d %d %d\n", c1->r, c1->g, c1->b, c2->r, c2->g, c2->b);
    return labDist(c1->r, c1->g, c1->b, c2->r, c2->g, c2->b);
}

typedef struct {
    unsigned char aveR;
    unsigned char aveG;
    unsigned char aveB;
    unsigned char medR;
    unsigned char medG;
    unsigned char medB;
    Colour* side;
    Colour** quantized;
    int count;
} ColourSummary;

typedef struct {
    int primaryIndex;
    int secondaryIndex;
    double dist;
    int minCount;
} Swap;


int comp_swap(const void* v1, const void* v2) {
    Swap x1=*((Swap*)v1);
    Swap x2=*((Swap*)v2);

    // return x1.dist - x2.dist;

    if(x1.dist < x2.dist) {
        return -1;
    } else if(x1.dist > x2.dist) {
        return 1;
    }

    // if(x1.minCount < x2.minCount) {
    //     return -1;
    // } else if(x1.minCount > x2.minCount) {
    //     return 1;
    // }

    return 0;
}

int comp_colour_count(const void* v1, const void* v2) {
    Colour* x1=*((Colour**)v1);
    Colour* x2=*((Colour**)v2);

    if(x1->count < x2->count) {
        return -1;
    } else if(x1->count > x2->count) {
        return 1;
    }

    return 0;
}

int comp_colour_metric_diff1(const void* v1, const void* v2) {
    Colour* x1=*((Colour**)v1);
    Colour* x2=*((Colour**)v2);

    if(x1->metricDiff < x2->metricDiff) {
        return -1;
    } else if(x1->metricDiff > x2->metricDiff) {
        return 1;
    } 

    return 0;
}

// int comp_colour_metric_diff2(const void* v1, const void* v2) {
//     Colour* x1=*((Colour**)v1);
//     Colour* x2=*((Colour**)v2);

//     if(x1->metricDiff2 < x2->metricDiff2) {
//         return -1;
//     } else if(x1->metricDiff2 > x2->metricDiff2) {
//         return 1;
//     } 

//     return 0;
// }

int comp_colour_binSearch(Colour* x1, Colour x2) {
    int cmp;
    if(x1->r < x2.r) {
        cmp = -1;
    } else if(x1->r > x2.r) {
        cmp = 1;
    } else {
        //==
        if(x1->g < x2.g) {
            cmp = -1;
        } else if(x1->g > x2.g) {
            cmp = 1;
        } else {
            //==
            if(x1->b < x2.b) {
                cmp = -1;
            } else if(x1->b > x2.b) {
                cmp = 1;
            } else {
                //==
                cmp = 0;
            }
        }
    }
    
    return cmp;
}

int comp_colour2(const void* v1, const void* v2) {
    Colour* x1=*((Colour**)v1);
    Colour* x2=*((Colour**)v2);

    double d1 = labDist(0, 0, 0, x1->r, x1->g, x1->b);
    double d2 = labDist(0, 0, 0, x2->r, x2->g, x2->b);

    if(d1 < d2) return -1;
    if(d1 > d2) return 1;

    return 0;
}

int comp_colour(const void* v1, const void* v2) {
    Colour* x1=*((Colour**)v1);
    Colour* x2=*((Colour**)v2);

    unsigned char c1[3];
    c1[0] = x1->r;
    c1[1] = x1->g;
    c1[2] = x1->b;

    unsigned char c2[3];
    c2[0] = x2->r;
    c2[1] = x2->g;
    c2[2] = x2->b;

    double* xyz1 = RGBtoXYZ(c1);
    double* xyz2 = RGBtoXYZ(c2);
    double* lab1 = XYZtoCIELab(xyz1);
    double* lab2 = XYZtoCIELab(xyz2);

    int cmp = 0;
    
    // if(false && lab1[0] < lab2[0]) {
    //     cmp = -1;
    // } else if(false && lab1[0] > lab2[0]) {
    //     cmp = 1;
    // } else {
    //     //==
    //     if(lab1[1] < lab2[1]) {
    //         cmp = -1;
    //     } else if(lab1[1] > lab2[1]) {
    //         cmp = 1;
    //     } else {
    //         //==
    //         if(lab1[2] < lab2[2]) {
    //             cmp = -1;
    //         } else if(lab1[2] > lab2[2]) {
    //             cmp = 1;
    //         } else {
    //             //==
    //             cmp = 0;
    //         }
    //     }
    // }

    if(cmp == 0 && lab1[0] < lab2[0]) cmp = -1;
    if(cmp == 0 && lab1[0] > lab2[0]) cmp = 1;
    if(cmp == 0 && lab1[1] < lab2[1]) cmp = -1;
    if(cmp == 0 && lab1[1] > lab2[1]) cmp = 1;
    if(cmp == 0 && lab1[2] < lab2[2]) cmp = -1;
    if(cmp == 0 && lab1[2] > lab2[2]) cmp = 1;
    
    free(xyz1);
    free(xyz2);
    free(lab1);
    free(lab2);

    return cmp;
}

int findTrueIndex(Colour** quantized, int index) {
    if(quantized[index]->trueIndex > 0 && quantized[index]->trueIndex != index) {
        return findTrueIndex(quantized, quantized[index]->trueIndex);
    }

    return index;
}

bool col_same(Colour* c1, Colour* c2) {
    return c1->r == c2->r && c1->g == c2->g && c1->b == c2->b;
}

int closest_colour(Colour** a, int n, unsigned char r, unsigned char g, 
        unsigned char b) {
    double smallest = INFINITY;
    int smallestKey = -1;
    double d = 0;
    for(int i = 0; i < n; i++) {
        d = labDist(r, g, b, a[i]->r, a[i]->g, a[i]->b);
        if(d < smallest) {
            smallest = d;
            smallestKey = i;
        }
    }

    return smallestKey;
}


void recalc_quantized_colours_counts(Colour** quantized, int num, int w, int h, int channels, unsigned char* image) {
    unsigned char* pixelOffset;
    unsigned char r, g, b;
    
    int found;

    for(int i = 0 ; i < num; i++) quantized[i]->count = 0;

    for(int i = 0; i < w; i++) {
        for(int j = 0; j < h; j++) {
            pixelOffset = image + (i + h * j) * channels;
            r = pixelOffset[0];
            g = pixelOffset[1];
            b = pixelOffset[2];

            found = closest_colour(quantized, num, r, g, b);
            if(found > -1)  {
                quantized[found]->count = quantized[found]->count + 1;
            } else {
                fputs("FAILUREXX", stderr);
                exit(1000);
            }
        }
    }
}

ColourSummary* calculateSummary(unsigned char* image, int w, int h, int channels, Colour** colours, int num) {
    if(num == 0) return NULL;
    printf("got %d distinct indexes\n", num);
    int desiredNum = num/20;
    if(desiredNum < 10) desiredNum = 10;
    if(desiredNum > 20) desiredNum = 20;
    // int desiredNum = 10;
    const int desiredSwaps = num - desiredNum;
    qsort(colours, num, sizeof(Colour*), comp_colour);

    int medNum = ceil((num + 1)/2);
    int countSoFar = 0;
    Colour* median = calloc(1, sizeof(Colour));

    long long int totalR;
    long long int totalG;
    long long int totalB;

    for(int i = 0; i < num; i++) {
        totalR += colours[i]->r;
        totalG += colours[i]->g;
        totalB += colours[i]->b;
        countSoFar += colours[i]->count;

        if(countSoFar <= medNum || (countSoFar - colours[i]->count + 1 
                <= medNum)) {
            memcpy(median, colours[i], sizeof(Colour));
        }
    }

    Colour* ave = calloc(1, sizeof(Colour));
    ave->r = totalR/num;
    ave->g = totalG/num;
    ave->b = totalB/num;

    Colour* distMetric = ave;

    double largestdiff = 1;

    for(int i = 0; i < num; i++) {
        if(colours[i]->count > 0) {
            double d = labDist(colours[i]->r, colours[i]->g, 
                colours[i]->b, distMetric->r, distMetric->g, distMetric->b);
            if(d > largestdiff) largestdiff = d;
        }
    }

    recalc_quantized_colours_counts(colours, num, w, h, channels, image);

    for(int i = 0; i < num; i++) {
        colours[i]->metricDiff = colours[i]->count * (1.5-labDist(colours[i]->r, colours[i]->g, 
                colours[i]->b, distMetric->r, distMetric->g, distMetric->b)/largestdiff);
    }


    Colour* side = calloc(1, sizeof(Colour));
    qsort(colours, num, sizeof(Colour*), comp_colour_metric_diff1);
    memcpy(median, colours[num - 1], sizeof(Colour));

    for(int i = 0; i < num; i++) {
        colours[i]->metricDiff *= labDist(colours[i]->r, colours[i]->g, 
                colours[i]->b, median->r, median->g, median->b)/largestdiff;
    }
    qsort(colours, num, sizeof(Colour*), comp_colour_metric_diff1);
    memcpy(side, colours[num - 1], sizeof(Colour));

    /** for debug**/
    qsort(colours, num, sizeof(Colour*), comp_colour);
    FILE* qFile = fopen("c_quantized.php", "w");
    fprintf(qFile, "<?php return [\n");
    for(int i = 0; i < num; i++) {
        fprintf(qFile, "['col' => [%d, %d, %d], 'n' => %d],\n", colours[i]->r, colours[i]->g, colours[i]->b, colours[i]->count);
    }
    fprintf(qFile, "]; ?>");
    fclose(qFile);
    /** for debug**/

    ColourSummary* ret = calloc(1, sizeof(ColourSummary));


    //DO SOMETHING WIHT Temp
    ret->aveR = ave->r;
    ret->aveG = ave->g;
    ret->aveB = ave->b;
    // ret->quantized = temp;
    ret->medR = median->r;
    ret->medG = median->g;
    ret->medB = median->b;
    ret->side = side;
    free(median);
    free(ave);
    return ret;
}

int search_colour(Colour** a, int n, unsigned char r, unsigned char g, 
        unsigned char b) {
    for(int i = 0; i < n; i++) {
        if(a[i]->r == r && a[i]->g == g && a[i]->b == b) return i;
    }

    return -1;
}

size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t written = fwrite(ptr, size, nmemb, stream);
    return written;
}

int main(int argc, char** argv) {
    hColoursMutex = CreateMutexW(NULL, FALSE, NULL);

    int x,y,channelCount;
    char* filename = argv[1];
    if(strlen(filename) == 0) {
        fprintf(stderr, "Usage: fmi filename\n");
        return 1;
    }

    if(strncmp(filename, "http", 4) == 0) {
        return 0;
        char* tempFile = "temp.png";

        // get_file(filename, tempFile);

        printf("using temp file - dl from internet\n");

        filename = tempFile;
    }

    unsigned char* data = stbi_load(filename, &x, &y, &channelCount, 0);
    unsigned char* data2 = calloc(MAX_SIZE*channelCount, sizeof(unsigned char));

    if(x <= 0 || y <= 0) {
        fprintf(stderr, "bad file");
        exit(1000);
    }

    if((x + y) > 2 * sqrt(MAX_SIZE)) {
        printf("orig %d x %d\n", x, y);
        double ratio = sqrt(MAX_SIZE)/(sqrt(x*y));
        int scaledX = x*ratio;
        int scaledY = y*ratio;
        printf("resizing image to %d x %d\n", scaledX, scaledY);
        if(!stbir_resize_uint8(data , x , y , 0,
                               data2, scaledX, scaledY, 0, channelCount)) {
            fprintf(stderr, "could not resize!\n");
            exit(10212);
        }
        x = scaledX;
        y = scaledY;
    }

    unsigned char* tmp = data;
    data = data2;
    free(tmp);

    if(data == NULL) {
        printf("could not read the image\n");
        return 1;
    }

    Colour** colours = calloc(x*y, sizeof(Colour*));

    int colIndex = 0;
    int dupes = 0;
    for(int i = 0; i < x; i++) {
        for(int j = 0; j < y; j++) {
            unsigned bytePerPixel = channelCount;
            unsigned char* pixelOffset = data + (i + y * j) * bytePerPixel;
            unsigned char r = pixelOffset[0];
            unsigned char g = pixelOffset[1];
            unsigned char b = pixelOffset[2];

            int found = search_colour(colours, colIndex, r, g, b);
            if(found > -1)  {
                dupes++;
                colours[found]->count++;
            } else {
                colours[colIndex] = calloc(1, sizeof(Colour));
                colours[colIndex]->r = r;
                colours[colIndex]->g = g;
                colours[colIndex]->b = b;
                colours[colIndex]->count = 1;
                colIndex++;
            }
            // fprintf(f, "%d => [%u, %u, %u], ", j, r, g, b);
        }
        // fprintf(f, "],\n");
    }

    // fprintf(f, "]; ?>");
    // fclose(f);
    // exit(1);

    printf("%d dupes %d total\n", dupes, x*y);

    ColourSummary* sum = calculateSummary(data, x, y, channelCount, colours, colIndex);
    if(sum == NULL) {
        printf("sum null\n");
        return 12;
    }
    printf("ave %d %d %d\n", sum->aveR, sum->aveG, sum->aveB);
    printf("med %d %d %d\n", sum->medR, sum->medG, sum->medB);

    printf("main %d %d %d\n", sum->medR, sum->medG, sum->medB);
    printf("side %d %d %d\n", sum->side->r, sum->side->g, sum->side->b);

    stbi_image_free(data);
}