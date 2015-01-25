#include"gucompute.h"

auto loadf(const char* fname, std::string pfx , bool bounded, int restarts){

    std::vector<gurunner> gs;

    std::ifstream fs(fname);
    if (!fs.good()) {
        std::cerr << " Input file "<< fname <<" is not good! "<<std::endl;
        exit(0);
    }
    std::string s;
    if (!std::getline(fs,s)) {
        std::cerr << "File Empty: "<<fname  <<std::endl;
    }
    while (std::getline(fs,s)){
        gurunner g = {pfx,"(PlaceHolder)", restarts, bounded, 0.5, 0.0, 1.0, 0.0, 1.0} ;
        g.load(s);
        g.train();
        gs.push_back(g);
    }
    return gs;
}

int main(int argc, char** argv){

    std::string prefix = "Test";
    unsigned num_threads=10;
    int restarts = 50;
    char* fname = "test_in.txt";

    //std::cerr << "(this program) fname prefix num_threads bounded num_retarts" << std::endl;
    //std::cerr << *argv[1] <<"," << *argv[2]<<std::endl;
    if (argc > 1) fname = argv[1];
    if (argc > 2) prefix = argv[2];
    //if (argc > 3) num_threads = (unsigned) *argv[3];
    //if (argc > 4) bounded_uniform = (bool)(int) *argv[4];
    //if (argc > 5) restarts = (int) *argv[5];
    //std::cerr << "Bound status: "<<bounded_uniform <<std::endl;
    //simple interface

    loadf(fname, prefix+"_bounded" , 1, restarts);
    loadf(fname, prefix+"_unbounded" , 0, restarts);
/*
 *    if (num_threads > gs.size()) num_threads = gs.size();
 *    std::vector <std::thread> thrs;
 *    thrs.resize(num_threads);
 *
 *    int idx = 0, idy = 0; // idx for counting runner objects idy for counting threads
 *    while (idx < gs.size()){
 *        if (idy<num_threads && idx < gs.size()) {
 *            auto g = gs[idx]; //current one to fork
 *            thrs[idy] = g.spawn();
 *            ++idy;
 *            ++idx;
 *        } else {
 *            //join after total number of threads spawn
 *            for (auto& th : thrs){
 *                th.join();
 *            }
 *            idy = 0;
 *        }
 *    }
 *    for (auto& th : thrs){
 *        th.join();
 *    } //finish up joining after jumping out
 */

    //for (auto & g: gs) g.train();

    return 0;
};
