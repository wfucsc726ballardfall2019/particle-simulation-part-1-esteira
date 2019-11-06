#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "common.h"
#define NEIGHBORS_SIZE 9

void getNeighbors(int pos, int num_bins, int *neighbors);
void clearNeighbors(int *neighbors);

using namespace std;

// Function to return which neighbors are valid
void clearNeighbors(int *neighbors) {
    for (int i = 0; i < NEIGHBORS_SIZE; i++) {
        neighbors[i] = -1;
    }
}

void getNeighbors(int pos, int num_bins, int *neighbors) {
    neighbors[0] = pos;

    if (pos % num_bins == 0) {
        // Left column
        if (floor(pos / num_bins) == 0)
        {
            // Top row => top left corner
            neighbors[1] = (pos + 1);
            neighbors[2] = (pos + num_bins);
            neighbors[3] = (pos + num_bins + 1);
            neighbors[4] = (-1);
            neighbors[5] = (-1);
            neighbors[6] = (-1);
            neighbors[7] = (-1);
            neighbors[8] = (-1);
        }
        else if (floor(pos / num_bins) == num_bins - 1)
        {
            // Bottom row => bottom left corner
            neighbors[1] = (pos + 1);
            neighbors[2] = (pos - num_bins);
            neighbors[3] = (pos - num_bins + 1);
            neighbors[4] = (-1);
            neighbors[5] = (-1);
            neighbors[6] = (-1);
            neighbors[7] = (-1);
            neighbors[8] = (-1);
        }
        else
        {
            neighbors[1] = (pos + 1);
            neighbors[2] = (pos - num_bins);
            neighbors[3] = (pos + num_bins);  
            neighbors[4] = (pos - num_bins + 1);
            neighbors[5] = (pos + num_bins + 1);    
            neighbors[6] = (-1);
            neighbors[7] = (-1);
            neighbors[8] = (-1);
        }
    }
    else if (pos % num_bins == num_bins - 1) {
        // Right column

        if (floor(pos / num_bins) == 0)
        {
            // Top row => top right corner
            neighbors[1] = (pos - 1);
            neighbors[2] = (pos + num_bins);
            neighbors[3] = (pos + num_bins - 1);
            neighbors[4] = (-1);
            neighbors[5] = (-1);
            neighbors[6] = (-1);
            neighbors[7] = (-1);
            neighbors[8] = (-1);
        }
        else if (floor(pos / num_bins) == num_bins - 1)
        {
            // Bottom row => bottom right corner
            neighbors[1] = (pos - 1);
            neighbors[2] = (pos - num_bins);
            neighbors[3] = (pos - num_bins - 1);
            neighbors[4] = (-1);
            neighbors[5] = (-1);
            neighbors[6] = (-1);
            neighbors[7] = (-1);
            neighbors[8] = (-1);
        }
        else 
        {
            neighbors[1] = (pos - 1);
            neighbors[2] = (pos - num_bins);
            neighbors[3] = (pos - num_bins - 1);
            neighbors[4] = (pos + num_bins);
            neighbors[5] = (pos + num_bins - 1);
            neighbors[6] = (-1);
            neighbors[7] = (-1);
            neighbors[8] = (-1);
        }
    }
    else if (floor(pos / num_bins) == 0)
    {
        // Top row
        neighbors[1] = (pos - 1);
        neighbors[2] = (pos + 1);
        neighbors[3] = (pos + num_bins);
        neighbors[4] = (pos + num_bins - 1);
        neighbors[5] = (pos + num_bins + 1);
        neighbors[6] = (-1);
        neighbors[7] = (-1);
        neighbors[8] = (-1);

    }
    else if(floor(pos / num_bins) == num_bins - 1)
    {
        // Bottom row
        neighbors[1] = (pos + 1);
        neighbors[2] = (pos - 1);
        neighbors[3] = (pos - num_bins);
        neighbors[4] = (pos - num_bins - 1);
        neighbors[5] = (pos - num_bins + 1);
        neighbors[6] = (-1);
        neighbors[7] = (-1);
        neighbors[8] = (-1);
    }
    else 
    {
        // All eight neighbors are valid
        neighbors[1] = (pos + 1);
        neighbors[2] = (pos - 1);
        neighbors[3] = (pos + num_bins);
        neighbors[4] = (pos - num_bins);
        neighbors[5] = (pos + num_bins + 1);
        neighbors[6] = (pos + num_bins - 1);
        neighbors[7] = (pos - num_bins + 1);
        neighbors[8] = (pos - num_bins - 1);
    }
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    // Set the size of the bins
    // The grid size is sqrt(0.0005 * number of particles)
    // So make sure the bins only consider the cutoff radius where the particles actually react to each other
    double bin_length = getCutoff() * 2;
    int num_bins = ceil((sqrt(getDensity() * n)) / bin_length);

    // Array to keep track of which particles are in which bins
    // We have to compute the bin for each particle, so this is O(n)
    vector<vector<particle_t> > bins(num_bins * num_bins);
    int offset_x;
    int offset_y;
    int which_bin;

    for (int i = 0; i < n; i++)
    {
        // Compute which bin a particle belongs to based on its location
        offset_x = floor(particles[i].x / bin_length);
        offset_y = floor(particles[i].y / bin_length);

        which_bin = num_bins * offset_y + offset_x;

        // Add the particle to the list of particles in that bin
        bins[which_bin].push_back(particles[i]);
    }

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    int *neighbors = new int[NEIGHBORS_SIZE];
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	    navg = 0;
        davg = 0.0;
	    dmin = 1.0;
        //
        //  compute forces
        //

        // Consider each particle
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;

            // Get the bin of the current particle
            offset_x = floor(particles[i].x / bin_length);
            offset_y = floor(particles[i].y / bin_length);

            which_bin = num_bins * offset_y + offset_x;

            // Get the neighbors of this bin 
            getNeighbors(which_bin, num_bins, neighbors);
            
            // Now we have valid neighbors, so compute force between the current particle and the particles in the neighboring bins
            // Consider each neighboring bin
            for (int k = 0; k < NEIGHBORS_SIZE; k++)
            {
                if (neighbors[k] >= 0) 
                {
                    // Consider each particle in that bin
                    for (int p = 0; p < bins[neighbors[k]].size(); p++)
                    {
                        // Compute the force between the current particle and the particles in this bin
                        apply_force(particles[i], bins[neighbors[k]][p], &dmin, &davg, &navg);
                    }
                }

            }

            // Clear the neighbors vector
            clearNeighbors(neighbors);
        }
 
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }

        // The particles have moved, so update the particles in each bin
        // First, clear the current bin information
        for (int i = 0; i < num_bins * num_bins; i++) 
            bins[i].resize(0);

        // Update
        for (int i = 0; i < n; i++)
        {
            // Compute which bin a particle belongs to based on its location
            offset_x = floor(particles[i].x / bin_length);
            offset_y = floor(particles[i].y / bin_length);

            which_bin = num_bins * offset_y + offset_x;

            // Add the particle to the list of particles in that bin

            bins[which_bin].push_back(particles[i]);
        }

        // Clear the neighbors vector
        clearNeighbors(neighbors);
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}