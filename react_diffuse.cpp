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

#define DefaultBin 1 // This represents only one side of the total bin
#define DefaultNumofRun 1
#define DefaultNumofParticleA 5  // initial population of a particles
#define DefaultNumofParticleB 5  // initial population of b particles
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

    list<int> a_particles;
    list<int> b_particles;

    int particlelocation_a[NUMofParticleA];
    int particlelocation_b[NUMofParticleB];
    int i, j, nstep;

    double reactiontime(double);

    ////////////////////////////////////////////////////
    // Parameter set
    ////////////////////////////////////////////////////

    double L = 1; // 1D interval length (half)

    double D = 1; // Diffusion rate

    double h = L / (BINNUM);    // interval length for the discretization
    double d = 2 * D / (h * h); // d is the jumping rate, including jumping to left or right
    double r = 1;               // reaction rate

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
    cout << "Diffusion rate = " << D << endl;
    cout << "jumping rate = " << d << endl;
    cout << "total iterations " << NUMofRUNS << endl;

    double exittime[NUMofRUNS];
    double total_reaction_time = 0.0;
    double total_jump_time = 0.0;
    int reaction_count = 0;
    int jump_count = 0;
    int c_population = 0; // Population of c particles

    for (int real = 0; real < NUMofRUNS; real++)
    {
        a_particles.clear(); // Clear previous run's data
        b_particles.clear();
        c_population = 0; 

        // Initialize a particles
        for (int i = 0; i < NUMofParticleA; i++)
        {
            a_particles.push_back(i); 
            particlelocation_a[i] = -L / 2; // a particles start at -L/2
        }

        // Initialize b particles
        for (int i = 0; i < NUMofParticleB; i++)
        {
            b_particles.push_back(i); 
            particlelocation_b[i] = L / 2; // b particles start at L/2
        }

        nstep = 0; 
        double timeTracker = 0.0; 
        bool exitflag = false; 
        cout << "Initial population of a particles: " << a_particles.size() << endl;
        cout << "Initial population of b particles: " << b_particles.size() << endl;
        cout << "Initial population of c particles: " << c_population << endl;
     
        while (!exitflag && (!a_particles.empty() || !b_particles.empty()))
        {   // total propensity = diffusion propensity + reaction propensity
            a0 = (d * (a_particles.size() + b_particles.size())) + (r * a_particles.size() * b_particles.size()); 

            tau = reactiontime(a0);
            timeTracker += tau;
            nstep++;

            r2 = 1.0 * rand() / RAND_MAX;
            double r2a0 = r2 * a0;

            if (r2a0 < d * a_particles.size())// a diffuse
            {
               
                int index = int(r2 * a_particles.size());
                auto it = a_particles.begin();
                advance(it, index);

                r2residual = r2 * a_particles.size() - index;
                if (r2residual > 0.5)
                    jumpdirection = 1;
                else
                    jumpdirection = -1;

                particlelocation_a[*it] += jumpdirection;
                if (particlelocation_a[*it] == BINNUM * jumpdirection)
                {
                    jumptimesfile << timeTracker << endl; 
                    total_jump_time += timeTracker;
                    a_particles.erase(it);
                }
            }
            else if (r2a0 < d * (a_particles.size() + b_particles.size()))// b diffuse
            {
                int index = int(r2 * b_particles.size());
                auto it = b_particles.begin();
                advance(it, index);

                r2residual = r2 * b_particles.size() - index;
                if (r2residual > 0.5)
                    jumpdirection = 1;
                else
                    jumpdirection = -1;

                particlelocation_b[*it] += jumpdirection;
                if (particlelocation_b[*it] == BINNUM * jumpdirection)
                {
                    jumptimesfile << timeTracker << endl; // Record jump time
                    total_jump_time += timeTracker;
                    b_particles.erase(it);
                }
            }
            else //reaction: a + b -> c
            {
                int a_index = int(r2 * a_particles.size());
                auto a_it = a_particles.begin();
                advance(a_it, a_index);

                int b_index = int(r2 * b_particles.size());
                auto b_it = b_particles.begin();
                advance(b_it, b_index);

                reactiontimefile << timeTracker << endl;
                total_reaction_time += timeTracker;
                reaction_count++;
                a_particles.erase(a_it);
                b_particles.erase(b_it);
                c_population++; // Increment c population
            }

            if (a_particles.empty() && b_particles.empty())
            {
                exitflag = true;
                exittime[real] = timeTracker;
                exittimefile << timeTracker << endl;
                exitstepfile << nstep << endl;
            }
        }

    }

    double sum = 0;
    for (i = 0; i < NUMofRUNS; i++)
        sum += exittime[i];
    double mean = sum / NUMofRUNS;
    cout << "final exiting time = " << sum / NUMofRUNS << endl;
    
    cout << "average reaction time " << total_reaction_time / reaction_count << endl;
    double variance = 0;
    for (i = 0; i < NUMofRUNS; i++)
        variance += (exittime[i] - mean) * (exittime[i] - mean);
    variance /= NUMofRUNS;


    cout << "end of simulation ..." << endl;
    cout << "final population of a particles: " << a_particles.size() << endl;
    cout << "final population of b particles: " << b_particles.size() << endl;
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