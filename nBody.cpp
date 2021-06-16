#include "omp.h"
#include <bits/types/FILE.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <unistd.h>
#include <chrono>

const double G_CONST = 6.67E-11;
const double DT = 10;

int nPart;

typedef struct {
    // X and Y positions
    double x, y;
    // Velocity
    double vx, vy;
    // Mass
    double m;
    // Force
    double fx, fy;
    // Acceleration
    //double a;
} Particles;


void InitParticles(int nPart, std::vector<Particles> &particles){
    std::srand(int(time(NULL)) * 1000);
    #pragma omp for
    for (auto &p : particles){
        p.m = drand48() * 1000 + 100;  // Mass
        p.vx = 0;
        p.vy = 0;
        p.x = drand48();  // x-position
        p.y = drand48();  // y-position
        p.fx = 0.0, p.fy = 0.0;
    }
}

double calcForces(Particles &i, Particles &j, bool x){
    double r_x = i.x - j.x;
    double r_y = i.y - j.y;
    double r = sqrt(r_x*r_x + r_y*r_y);

    if (x == true){
        if (r_x == 0){
            return 0;
        }
        return (G_CONST * i.m * j.m) / (r*r*r) * r_x;
    }
    else {
        if (r_y == 0){
            return 0;
        }
        return (G_CONST * i.m * j.m) / (r*r*r) * r_y;
    }
}

int main(int argc, char *argv[]){
    int t_step = atoi(argv[2]);
    nPart = atoi(argv[1]);
    std::vector<Particles> particles(nPart, Particles()); // Particles on all nodes
    int totpart;
    double sim_t, startTime, endTime;
    FILE *fp, *fp2;
    char wfile[1024];
    sprintf(wfile, "timedat.%d",0);
    fp=fopen(wfile, "w+");

    omp_lock_t lock;
    omp_init_lock(&lock);

    auto start = std::chrono::steady_clock::now();

    InitParticles(nPart, particles);

    #pragma omp parallel
    for (int i=0; i<nPart; i++){
        fprintf(fp,"%d %d %d %f %f %f\n", 0, i, t_step, particles[i].x, particles[i].y, particles[i].m);
	    printf("%d %d %f %f %f\n", i, 0, particles[i].x, particles[i].y, particles.at(i).m);
	}
	fprintf(fp,"\n\n");

    
    for (int i = 1; i < t_step; i++){

        #pragma omp for
        for (int l = 0; l < particles.size(); l++){
            particles.at(l).fx = 0.0;
            particles.at(l).fy = 0.0;
        }
       
        #pragma omp for
        for (int m = 0; m < nPart; m++){
            for (int k = m+1; k < nPart; k++){
                double fx_1 = calcForces(particles.at(m), particles.at(k), true);
                double fy_1 = calcForces(particles.at(m), particles.at(k), false); 
                omp_set_lock(&lock);
                particles.at(m).fx -= fx_1;
                particles.at(m).fy -= fy_1; 
                omp_unset_lock(&lock);
            }                
        }

        #pragma omp for
        for (Particles &p : particles){
            p.x += DT*p.vx;
            p.y += DT*p.vy;
            p.vx += DT/p.m*p.fx;
            p.vy += DT/p.m*p.fy;
        }

        #pragma omp single
        for (int j=0; j<nPart; j++){
            fprintf(fp,"%d %d %d %f %f %f\n",i, j, t_step, particles[j].x, particles[j].y, particles[j].m);
            printf("%d %d %f %f %f\n", i, j, particles[j].x, particles[j].y, particles.at(j).m);
        }
        fprintf(fp, "\n\n");
    }

    auto end = std::chrono::steady_clock::now();

    
    std::cout << "Elapsed time in milliseconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;
 
    std::cout << "Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        << " sec";

    fclose(fp);
    omp_destroy_lock(&lock);

    
    return 0;
}