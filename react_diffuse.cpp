#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <random>
#include <list>
#include <iterator>

using namespace std;

#define DefaultBin 4// This represents only one side of the total bin
#define DefaultNumofRun 100
#define DefaultNumofParticleA 1 // initial population of a particles
#define DefaultNumofParticleB 1 // initial population of b particles
// #define printtraj
// #define DEBUG

int main(int argc, char *argv[])
{
    int BINNUM, NUMofRUNS, NUMofParticleA, NUMofParticleB;
    if (argc > 1)
        BINNUM = atoi(argv[1]);
    else
        BINNUM = DefaultBin;
    if (BINNUM == 0)
        BINNUM = DefaultBin;

    if (argc > 2)
        NUMofRUNS = atoi(argv[2]);
    else
        NUMofRUNS = DefaultNumofRun;
    if (NUMofRUNS == 0)
        NUMofRUNS = DefaultNumofRun;

    if (argc > 3)
        NUMofParticleA = atoi(argv[3]);
    else
        NUMofParticleA = DefaultNumofParticleA;
    if (NUMofParticleA == 0)
        NUMofParticleA = DefaultNumofParticleA;

    if (argc > 4)
        NUMofParticleB = atoi(argv[4]);
    else
        NUMofParticleB = DefaultNumofParticleB;
    if (NUMofParticleB == 0)
        NUMofParticleB = DefaultNumofParticleB;


    double particlelocation_a[NUMofParticleA];
    double particlelocation_b[NUMofParticleB];
    int i, j, nstep;

    double reactiontime(double);

    ////////////////////////////////////////////////////
    // Parameter set
    ////////////////////////////////////////////////////

    double L = 1; // 1D interval length (half)

    double D = 1; // Diffusion rate

    double h = L / (BINNUM);    // interval length for the discretization
    double d = 2 * D / (h * h); // d is the jumping rate, including jumping to left or right
    double k = 1;
    double r = k / h; // reaction rate

    //*****************************************
    // variables for SSA process
    //*****************************************
    double a0;
    double r1, r2, r2residual;
    double tau;
    int jump_index, jumpdirection;

    srand(1.0);

    ofstream exittimefile("exittimefile", ios::out);
    ofstream reactiontimefile("reactiontimes", ios::out);
    ofstream jumptimesfile("jumptimes", ios::out);
    ofstream exitstepfile("exitstep", ios::out);

    //************************************************
    // Begin simulation, and set up timer
    //************************************************
    cout << "Begin the model simulation ..." << endl;
    cout << "BINNUM = " << BINNUM << endl;
    // cout << "Diffusion rate = " << D << endl;
    cout << "jumping rate = " << d << endl;
    cout << "reacting rate = " << r << endl;
    cout << "total iterations " << NUMofRUNS << endl;
    // NEW FILE to store times when a- and b-particles meet but do not react
    ofstream meetfile("meet_times.txt", ios::out);
    double exittime[NUMofRUNS];
    double total_reaction_time = 0.0;
    double total_jump_time = 0.0;
    int reaction_count = 0;
    int jump_count = 0;
    int c_population = 0; // Population of c particles
    int a_jump, b_jump;
    for (int real = 0; real < NUMofRUNS; real++)
    {
        c_population = 0;

        // Initialize a particles
        for (int i = 0; i < NUMofParticleA; i++)
        {
            particlelocation_a[i] = -BINNUM / 2; // a particles start at -BINNUM/2
        }

        // Initialize b particles
        for (int i = 0; i < NUMofParticleB; i++)
        {
            particlelocation_b[i] = BINNUM / 2; // b particles start at BINNUM/2

        }

        nstep = 0;
        double timeTracker = 0.0;
        bool exitflag = false;
        
        double reactionrate;
        int samelocation = 0;
     
        while (!exitflag)
        {  
            reactionrate = r*samelocation;
            // total propensity = diffusion propensity + reaction propensity
            a0 = d*(NUMofParticleA + NUMofParticleB) + reactionrate;
            tau = reactiontime(a0);

            timeTracker += tau;
            nstep++;

            r2 = 1.0 * rand() / RAND_MAX;
            double r2residual = r2 * a0;
            if (r2residual < d * (NUMofParticleA + NUMofParticleB)) // Diffusion
            {
                if (r2residual < d*NUMofParticleA)// A jump
                {    
                    a_jump+=1;  
                    r2residual /= d;

                    int index = int(r2residual);
                    r2residual   =r2residual-index ; 
    
                    if (r2residual > 0.5)
                        jumpdirection = 1;
                    else
                        jumpdirection = -1;
                    // Reflective boundary. If particle jump out of boundary, jump backward instead
                    if (particlelocation_a[index] == BINNUM * jumpdirection)
                    {
                        particlelocation_a[index] -= jumpdirection;
                    }
                    else
                    {
                        particlelocation_a[index] += jumpdirection;
                    }
                    for (int i=0; i < NUMofParticleB;i++){
                        
                    if (particlelocation_a[index] == particlelocation_b[i])
                    {
                        samelocation = 1;
                    }
                    else{
                        samelocation = 0;
                    }
     
                    }
                   
                }
                else
                {
    
                    b_jump+=1;
                    r2residual = (r2residual - d*NUMofParticleA)/d;

                    int index = int(r2residual);
                    r2residual   =r2residual-index; 

                    if (r2residual > 0.5)
                        jumpdirection = 1;
                    else
                        jumpdirection = -1;

                    // Reflective boundary. If particle jump out of boundary, jump backward instead
                    if (particlelocation_b[index] == BINNUM * jumpdirection)
                    {
                        particlelocation_b[index] -= jumpdirection;
                    }
                    else
                    {
                        particlelocation_b[index] += jumpdirection;
                    }
                    if (particlelocation_a[index] == particlelocation_b[index])
                    {
                        samelocation = 1;
                    }
                    else{
                        samelocation = 0;
                    }
                }
              
            
            }
            else // reaction: a + b -> c
            {   
            reaction_count +=1;
            //     printf("total_reaction_time before is %f\n",total_reaction_time);
            //     printf("timeTracker is %f\n",timeTracker);
             total_reaction_time += timeTracker;
            //     printf("total_reaction_time after is %f\n",total_reaction_time);
            //     int a_index = int(r2 * a_particles.size());
            //     auto a_it = a_particles.begin();
            //     advance(a_it, a_index);

            //     int b_index = int(r2 * b_particles.size());
            //     auto b_it = b_particles.begin();
            //     advance(b_it, b_index);
               
            //     a_particles.erase(a_it);
            //     b_particles.erase(b_it);
            //     c_population++;

            exitflag = true;
            }

          
        }
    }

    cout << "average reaction time " << total_reaction_time / reaction_count << endl; 
    printf("Number of jumps of A before reaction fire %d\n", a_jump);
    printf("Number of jumps of B before reaction fire %d\n", b_jump);
    cout << "end of simulation ..." << endl;
    cout << "final population of c particles: " << c_population << endl;
}


double reactiontime(double a0)
{
    double r1, tau;

    // reaction time follows exponential
    do
    {
        r1 = 1.0 * rand() / RAND_MAX;
    } while (r1 <= 0 || r1 >= 1); // random number
    tau = -1.0 / a0 * log(r1);

    return tau;
}

