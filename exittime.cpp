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

using namespace std;

#define DefaultBin 50 //This represent only one side of the total bin
#define DefaultNumofRun 1e6
#define DefaultNumofParticle 2
// #define printtraj
// #define DEBUG

int main(int argc, char* argv[]){

  int BINNUM, NUMofRUNS, NUMofParticle;
  if(argc>1) 
	BINNUM=atoi(argv[1]);
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
    
	int particlenumber = 1;

	int particlelocation[NUMofParticle];
    int i, j, nstep;

    double reactiontime(double );

  ////////////////////////////////////////////////////
  //Parameter set
  ////////////////////////////////////////////////////

    double L = 1; // 1D interval length (half) 

    double D = 1; // Diffusion rate
    
    double h = L/(BINNUM);  // interval length for the discretization
    double d = 2*D/(h*h); // d is the jumping rate, including jumping to left or right
   
  //*****************************************
  //variables for SSA process
  //*****************************************
  double a0;
    a0 = NUMofParticle * d; // the total jumping propensity, each particle has the same jumping rate, so the calculation is easy.
  double r1, r2, r2residual;
  double tau;
    int jump_index, jumpdirection;
 
    srand(time(NULL));
    //srand(1.0); 
  
    ofstream exittimefile("exittimefile", ios::out);
    ofstream exitstepfile("exitstep", ios::out);

  //************************************************
  //Begin simulation, and set up timer
  //************************************************
    cout<<"Begin the model simulation ..."<<endl;
    cout<<"BINNUM = "<< BINNUM << endl;
    cout << "initial number of particles: " << NUMofParticle << endl;
    cout<<"Diffusion rate = " << D <<endl;
    cout << "jumping rate = " << d << endl;
    cout<<"total iterations "<< NUMofRUNS <<endl;
	 double exittime[NUMofRUNS];

  for(int real=0; real< NUMofRUNS; real++){
      
      for (i=0; i < NUMofParticle; i++)
          particlelocation[i] = 0; // reset particle locations. all particles are at location 0
      nstep = 0;  // reset number of jumping steps
      
      double timeTracker = 0.0; // reset time
	    bool exitflag = false;  // reset exitflag

      while(!exitflag){
          
          tau = reactiontime(a0);
          
          timeTracker += tau;
          nstep++;

      //***********************************

		  r2 = 1.0*rand()/RAND_MAX;
		  jump_index = (int) (r2* NUMofParticle); // uniform distribution, select the index for the particle that will jump
          
          

      r2residual = r2*NUMofParticle - jump_index; //decide left or right jump
      if (r2residual > 0.5)
          jumpdirection = 1;
      else
          jumpdirection = -1;

          
      // cout << "one step" << r2residual << "jump" << jumpdirection << endl;
          
		  particlelocation[jump_index] += jumpdirection;
          
      // cout << jump_index << "jump to" << particlelocation[jump_index] << endl;
          
      if (particlelocation[jump_index] == BINNUM*jumpdirection) {
          exitflag = true;
          // cout << "time = " << timeTracker << endl;
          exittimefile << timeTracker << endl;
                exitstepfile << nstep << endl;
          exittime[real] = timeTracker;
          // cout << timeTracker << endl;
		  }
	  }
    }
    
	double sum = 0;
	for (i = 0; i < NUMofRUNS; i++) 
		sum += exittime[i];
	cout << "average exiting time = " << sum/NUMofRUNS << endl;
    cout<<"end of simulation ..."<<endl;
}

double reactiontime(double a0) {
    double r1, tau;

    // reaction time follows exponential
    // do {r1=1.0*rand()/RAND_MAX;} while(r1<=0 || r1>=1); //random number
    // tau = -1.0/a0*log(r1);

    // fixed value
    // double tau = 1/a0;  


    //uniform distribution
    //  do {r1=1.0*rand()/RAND_MAX;} while(r1<=0 || r1>=1);
    //  tau = 2/a0*r1;
       
    
    //try normal distribution here
    static std::default_random_engine generator;  
    static std::normal_distribution<double> distribution(1.0, 0.5); // mean = 1.0, stddev = 0.5
    do {
           r1 = distribution(generator); 
        } while (r1 <= 0); 
    tau = r1 / a0; //mean = 1.0/a0, stddev = 0.5/a0
    
    return tau;
}
