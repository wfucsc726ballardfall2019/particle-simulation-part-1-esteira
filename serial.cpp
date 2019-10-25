#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "common.h"

using namespace std;

// Function to return which neighbors are valid
vector<int> getNeighbors(int pos, int num_bins) {
    //cout << "** Considering pos " << pos << endl;
    vector<int> neighbors(0);
    if (pos % num_bins == 0) {
        //cout << "A" << endl;
        // Left column
        if (floor(pos / num_bins) == 0)
        {
            //cout << "AA" << endl;
            // Top row => top left corner
            neighbors.push_back(pos);
            neighbors.push_back(pos + 1);
            neighbors.push_back(pos + num_bins);
            neighbors.push_back(pos + num_bins + 1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
        }
        else if (floor(pos / num_bins) == num_bins - 1)
        {
            //cout << "AB" << endl;
            // Bottom row => bottom left corner
            neighbors.push_back(pos);
            neighbors.push_back(pos + 1);
            neighbors.push_back(pos - num_bins);
            neighbors.push_back(pos - num_bins + 1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
        }
        else
        {
            //cout << "AC" << endl;
            neighbors.push_back(pos);
            neighbors.push_back(pos + 1);
            neighbors.push_back(pos - num_bins);
            neighbors.push_back(pos + num_bins);  
            neighbors.push_back(pos - num_bins + 1);
            neighbors.push_back(pos + num_bins + 1);    
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
        }
    }
    else if (pos % num_bins == num_bins - 1) {
        // Right column
        //cout << "B" << endl;

        if (floor(pos / num_bins) == 0)
        {
            //cout << "BA" << endl;
            // Top row => top right corner
            neighbors.push_back(pos);
            neighbors.push_back(pos - 1);
            neighbors.push_back(pos + num_bins);
            neighbors.push_back(pos + num_bins - 1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
        }
        else if (floor(pos / num_bins) == num_bins - 1)
        {
            //cout << "BB" << endl;
            // Bottom row => bottom right corner
            neighbors.push_back(pos);
            neighbors.push_back(pos - 1);
            neighbors.push_back(pos - num_bins);
            neighbors.push_back(pos - num_bins - 1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
        }
        else 
        {
            //cout << "BC" << endl;
            neighbors.push_back(pos);
            neighbors.push_back(pos - 1);
            neighbors.push_back(pos - num_bins);
            neighbors.push_back(pos - num_bins - 1);
            neighbors.push_back(pos + num_bins);
            neighbors.push_back(pos + num_bins - 1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
            neighbors.push_back(-1);
        }
    }
    else if (floor(pos / num_bins) == 0)
    {
        //cout << "C" << endl;
        // Top row
        neighbors.push_back(pos);
        neighbors.push_back(pos - 1);
        neighbors.push_back(pos + 1);
        neighbors.push_back(pos + num_bins);
        neighbors.push_back(pos + num_bins - 1);
        neighbors.push_back(pos + num_bins + 1);
        neighbors.push_back(-1);
        neighbors.push_back(-1);
        neighbors.push_back(-1);

    }
    else if(floor(pos / num_bins) == num_bins - 1)
    {
        //cout << "D" << endl;
        // Bottom row
        neighbors.push_back(pos);
        neighbors.push_back(pos + 1);
        neighbors.push_back(pos - 1);
        neighbors.push_back(pos - num_bins);
        neighbors.push_back(pos - num_bins - 1);
        neighbors.push_back(pos - num_bins + 1);
        neighbors.push_back(-1);
        neighbors.push_back(-1);
        neighbors.push_back(-1);
    }
    else 
    {
        //cout << "E" << endl;
        // All eight neighbors are valid
        neighbors.push_back(pos);
        neighbors.push_back(pos + 1);
        neighbors.push_back(pos - 1);
        neighbors.push_back(pos + num_bins);
        neighbors.push_back(pos - num_bins);
        neighbors.push_back(pos + num_bins + 1);
        neighbors.push_back(pos + num_bins - 1);
        neighbors.push_back(pos - num_bins + 1);
        neighbors.push_back(pos - num_bins - 1);
    }
    /*
    cout << "**After making neighbors list size is " << neighbors.size() << endl;
    cout << "**Elements in there..." << endl;
    for (int i = 0; i < neighbors.size(); i++)
    {
        cout << neighbors[i] << " ";
    }
    cout << endl;
    */
    return neighbors;
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
    // The grid size is 0.0005 * number of particles
    // So make sure the bins only consider the cutoff radius where the particles actually react to each other
    // We know cutoff is the length of the bin
    //int num_bins = ceil((sqrt(getDensity() * n)) / getCutoff());
    double bin_length = getCutoff() * 2;
    int num_bins = ceil((sqrt(getDensity() * n)) / bin_length);
    cout << "num_bins " << num_bins << endl;

    // Array to keep track of which particles are in which bins
    // We have to compute the bin for each particle, so this is O(n)
    vector<vector<particle_t> > bins(num_bins * num_bins);
    int offset_x;
    int offset_y;
    int which_bin;

    for (int i = 0; i < n; i++)
    {
        // Compute which bin a particle belongs to based on its location
        //offset_x = floor(particles[i].x / getCutoff());
        //offset_y = floor(particles[i].y / getCutoff());
        offset_x = floor(particles[i].x / bin_length);
        offset_y = floor(particles[i].y / bin_length);

        which_bin = num_bins * offset_y + offset_x;

        //cout << "particle: " << i << ", xpos: " << particles[i].x << ", xoff: " << offset_x << ", ypos: " << particles[i].y << ",yoff: " << offset_y << ", bin: " << which_bin << endl;

        if (which_bin >= num_bins * num_bins) {
            cout << "out of boundaries" << endl;
            cout << "w " << which_bin << endl;
            cout << "x " << offset_x << endl;
            cout << "y " <<offset_y << endl;
            cout << "x pos: " << particles[i].x << endl;
            cout << "y pos: " << particles[i].y << endl;
        }

        // Add the particle to the list of particles in that bin
        bins[which_bin].push_back(particles[i]);
    }

    /*
    cout << "## ORIGINAL BINS ##" << endl;
    for (int t = 0; t < bins.size(); t++)
    {
        cout << "In bin " << t << endl;
        for (int s = 0; s < bins[t].size(); s++)
        {
            cout << "particle x:" << bins[t][s].x << " particle y:" << bins[t][s].y << endl;
        }
    }
    */


    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    vector<int> neighbors(9);
	
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
            //offset_x = floor(particles[i].x / getCutoff());
            //offset_y = floor(particles[i].y / getCutoff());
            offset_x = floor(particles[i].x / bin_length);
            offset_y = floor(particles[i].y / bin_length);

            which_bin = num_bins * offset_y + offset_x;

            // Get the neighbors of this bin 
            //cout << "particle " << i << " neighbors" << endl;
            neighbors = getNeighbors(which_bin, num_bins);
            /*
            for (int t = 0; t < neighbors.size(); t++){
                if (neighbors[t] >= 0){
                    cout << neighbors[t] << " ";
                }
            }
            cout << endl;
            */
            

            // Now we have valid neighbors, so compute force between the current particle and the particles in the neighboring bins
            // Consider each neighboring bin
            for (int k = 0; k < neighbors.size(); k++)
            {
                if (neighbors[k] >= 0) 
                {
                    /*
                    cout << "~~~Considering bin: " << neighbors[k] << endl;
                    cout << "~~~Number particles in that bin: " << bins[neighbors[k]].size() << endl;
                    cout << "~~~ BINS ~~~" << endl;
                    for (int t = 0; t < bins.size(); t++)
                    {
                        cout << "In bin " << t << endl;
                        for (int s = 0; s < bins[t].size(); s++)
                        {
                            cout << "particle x:" << bins[t][s].x << " particle y:" << bins[t][s].y << endl;
                        }
                    }
                    */
                    // Consider each particle in that bin
                    for (int p = 0; p < bins[neighbors[k]].size(); p++)
                    {
                        //cout << "~~~particle x " << bins[neighbors[k]][p].x << " particle y" << bins[neighbors[k]][p].y << endl;
                        // Compute the force between the current particle and the particles in this bin
                        if (bins[neighbors[k]].size() > 0)
                        {
                            apply_force(particles[i], bins[neighbors[k]][p], &dmin, &davg, &navg);
                        }
                    }
                }

            }

            // Clear the neighbors vector
            neighbors.clear();
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
            //offset_x = floor(particles[i].x / getCutoff());
            //offset_y = floor(particles[i].y / getCutoff());
            offset_x = floor(particles[i].x / bin_length);
            offset_y = floor(particles[i].y / bin_length);

            which_bin = num_bins * offset_y + offset_x;

            // Add the particle to the list of particles in that bin

            //cout << "This is the particle at x-axis and y-axis: " << particles[i].x << particles[i].y << endl;
            //cout << "Belonging to bin #" << which_bin << endl;
            bins[which_bin].push_back(particles[i]);
        }
        /*
        cout << "## AFTER UPDATING BINS ##" << endl;
        for (int t = 0; t < bins.size(); t++)
        {
            cout << "In bin " << t << endl;
            for (int s = 0; s < bins[t].size(); s++)
            {
                cout << "particle x:" << bins[t][s].x << " particle y:" << bins[t][s].y << endl;
            }
        }
        */

        // Clear the neighbors vector
        neighbors.clear();
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

// Original 
/*
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

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
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	    navg = 0;
        davg = 0.0;
	    dmin = 1.0;
        //
        //  compute forces
        //
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
				apply_force( particles[i], particles[j],&dmin,&davg,&navg);
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
*/
