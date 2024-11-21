//**********************************************
// This piece of code intend to test the exit time for one or multiple particles exit from a 1D domain [-L L]. Here we assume that the particles are located in origin. They can be modelled by jumps in discrete grids or a Gaussian diffusion equation.
//
//**********************************************

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <random>
#include <algorithm>
#include <list>
#include <iterator>

using namespace std;

#define DefaultBin 1// This represent only one side of the total bin
#define DefaultNumofRun 1
#define DefaultNumofParticle 10000
// #define printtraj
// #define DEBUG

int main(int argc, char *argv[])
{

  int BINNUM, NUMofRUNS, NUMofParticle;
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
    NUMofParticle = atoi(argv[3]);
  else
    NUMofParticle = DefaultNumofParticle;
  if (NUMofParticle == 0)
    NUMofParticle = DefaultNumofParticle;


  list<int> a_particles;
  list<int> b_particles; 

  int particlelocation[NUMofParticle];
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
  // a0 = (d + r)*a_particles.size(); // total propensityle * d; // the total jumping propensity, each particle has the same jumping rate, so the calculation is easy.
  //a0=d+r;
  double r1, r2, r2residual;
  double tau;
  int jump_index, jumpdirection;

  //srand(time(NULL));
  srand(1.0);

  ofstream exittimefile("exittimefile", ios::out);
  ofstream exitstepfile("exitstep", ios::out);

  //************************************************
  // Begin simulation, and set up timer
  //************************************************
  cout << "Begin the model simulation ..." << endl;
  cout << "BINNUM = " << BINNUM << endl;
  //cout << "initial number of particles: " << NUMofParticle << endl;
  cout << "Diffusion rate = " << D << endl;
  cout << "jumping rate = " << d << endl;
  cout << "total iterations " << NUMofRUNS << endl;
  double exittime[NUMofRUNS];
  for (int real = 0; real < NUMofRUNS; real++)
  {  
        
        a_particles.clear();//Clear previous run's data
        b_particles.clear();
        for (int i = 0; i < NUMofParticle; i++) {
          a_particles.push_back(i); //reset a particle
          particlelocation[i] = 0; // reset particle locations. all particles are at location 0
        }
        nstep = 0; // reset number of jumping steps
        double timeTracker = 0.0; // reset time
        bool exitflag = false; // reset exitflag
        cout << "Initial population of a particles: " << a_particles.size() << endl;
        cout << "Initial population of b particles: " << b_particles.size() << endl;
        while (!exitflag && !a_particles.empty()) {

            a0 = (d + r)*a_particles.size(); //the total jumping and reacting propensity

            tau = reactiontime(a0);
            timeTracker += tau;
            nstep++;

            r2 = 1.0 * rand() / RAND_MAX;
            double r2a0 = r2 * a0;

            int index = int(r2* a_particles.size());//select the index for the particle that will jump or react
            auto it = a_particles.begin();//initialize iterator at the begining of the a_particles list 
            advance(it, index); // Move iterator point to the selected index
   
            if (r2a0 < d) { // Diffusion event
                r2residual = r2 * a_particles.size() - index; // decide left or right jump
                if (r2residual > 0.5)
                    jumpdirection = 1;
                else
                    jumpdirection = -1;
               
                particlelocation[*it] += jumpdirection;
                if (particlelocation[*it] == BINNUM * jumpdirection) {
                    a_particles.erase(it); // Remove particle from a_particles list
                }
            } else { // Reaction event
               
              a_particles.erase(it);// Remove particle from a_particles list
              b_particles.push_back(*it); //Add particle to b_particles list
            } 

            if (a_particles.empty()) {
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
  cout << "average exiting time = " << sum / NUMofRUNS << endl;

  double variance = 0;
  for (i = 0; i < NUMofRUNS; i++)
    variance += (exittime[i] - mean) * (exittime[i] - mean);
  variance /= NUMofRUNS;
  //cout << "Variance of exiting time = " << variance << endl;
  cout << "end of simulation ..." << endl;
  cout << "final population of a particles: " << a_particles.size()   << endl;
  cout << "final population of b particles: " << b_particles.size()   << endl;
   
}

double reactiontime(double a0)
{
  double r1, tau;

  // reaction time follows exponential
  // do {r1=1.0*rand()/RAND_MAX;} while(r1<=0 || r1>=1); //random number
  // tau = -1.0/a0*log(r1);

  // fixed value
  // double tau = 1/a0;

  // uniform distribution
  //   do {r1=1.0*rand()/RAND_MAX;} while(r1<=0 || r1>=1);
  //   tau = 2/a0*r1;

  // try normal distribution here
  static std::default_random_engine generator;
  static std::normal_distribution<double> distribution(1.0, 0.5); // mean = 1.0, stddev = 0.5
  do
  {
    r1 = distribution(generator);
  } while (r1 <= 0);
  tau = r1 / a0; // mean = 1.0/a0, stddev = 0.5/a0

  return tau;
}