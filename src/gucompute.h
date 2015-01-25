/*
    Gaussian - Uniform mixture trainer written in c++11

    One runner is associated with a simple set of data points and variables.
    It is trained multiple times with random restarts and then the best values are
    written out.

    Yu Liu 2015
 */
#ifndef gucompute_h
#define gucompute_h

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<cmath>
#include<algorithm>
#include<random>
#define pi M_PI
#define LOG2 1
#if LOG2
#define LOG 0
#endif
static double pi2 = sqrt(2*pi) ;

class gurunner
{
    public:
        gurunner (const char* pfx, const char*g,  const unsigned& num,
                const bool& uu,const double& ww,const double& aa,const double& bb,
                const double& muu, const double& sigmaa)
            : prefix (pfx), gene(g),uniform_fixed(uu), restart_num(num), w(ww),
            a(aa), b(bb), mu(muu), sigma(sigmaa), res(-1.0), current_run_idx(0)
        {
            as.resize(restart_num);
            bs.resize(restart_num);
            ws.resize(restart_num);
            mus.resize(restart_num);
            sigmas.resize(restart_num);
            ress.resize(restart_num);
        };

        gurunner ()  : gurunner ("Test","TestGene", 5, 0, 0.5, 0.0, 1.0, 0.0, 1.0) {};
        //Delegation

        ~gurunner() {};

        void test ()
        {
            std::cout << "Testing \t "<<prefix << " " << pi2 <<std::endl;
        };

        void load (const char*);
        void init () ; // randomly initialize
        void record (); // record current parameters
        void writeout(const unsigned&); // write out a report of the whole run
        void run (double, int) ; //main routine
        void train (); //run with restarts


        //running routines
        void get_member_likelihood(); //generate uniform and normal likelihoods based on the weights and params
        void get_w(); //calculate new w;
        void get_uni_params(); //if uniform is not bounded at lower/higher datum, calculate the new bound.
        void get_normal_params(); // calculate the sigma and mu given the datum.
        void get_total_likelihood(); // calculate res

    private:
        std::string prefix, gene;
        bool uniform_fixed;
        //whether the uniform distribution is bounded lower at 0 always
        double a, b; //initial values for U(a,b) uniform distribution of the mixture
        double w; //weight prior
        double mu, sigma; // Gaussian parameters
        const unsigned restart_num; //number of random restarts
        unsigned sz_data, current_run_idx;
        double res; // Likelihood value
        std::vector<double> as, bs, ws, mus, sigmas, ress; // container for restarts
        std::vector<double> data, Pg, Pu, W;
        //constexpr double pi2=sqrt(pi); //valid in D not in C++11

        std::default_random_engine rndgen;

};



#endif
