#include"gucompute.h"

auto routine( std::string fname, std::string pfx , bool bounded, int restarts){

    std::vector<std::shared_ptr<gurunner>> gs;

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
        auto  g =  std::shared_ptr<gurunner>(new gurunner {pfx,"(PlaceHolder)", restarts, bounded, 0.5, 0.0, 1.0, 0.0, 1.0} );
        g->load(s);
        //g->train();
        gs.push_back(g);
    }
    return gs;
}

int main(int argc, char** argv){

    std::string prefix = "Test";
    unsigned num_threads=10, restarts = 500;
    std::string fname = "test_in.txt";

    //std::cerr << "(this program) fname prefix num_threads bounded num_retarts" << std::endl;
    //std::cerr << *argv[1] <<"," << *argv[2]<<std::endl;
    if (argc > 1) fname = argv[1];
    if (argc > 2) prefix = argv[2];
    //if (argc > 3) num_threads = (unsigned) *argv[3];
    //if (argc > 4) bounded_uniform = (bool)(int) *argv[4];
    //if (argc > 5) restarts = (int) *argv[5];
    //std::cerr << "Bound status: "<<bounded_uniform <<std::endl;
    //simple interface

    auto gs = routine(fname, prefix+"_bounded" , 1, restarts);
    routine(fname, prefix+"_unbounded" , 0, restarts);
    if (num_threads > gs.size()) num_threads = gs.size();
    std::vector <std::thread> thrs;
    thrs.resize(num_threads);

    unsigned idx = 0, idy = 0; // idx for counting runner objects idy for counting threads
    while (idx < gs.size()){
        if (idy<num_threads && idx < gs.size()) {
            auto g = gs[idx]; //current one to fork
            thrs[idy] = g->spawn();
            ++idy;
            ++idx;
        } else {
            //join after total number of threads spawn
            for (unsigned idz =0; idz<idy; ++idz){
                thrs[idz].join();
            }
            std::cout<< "Joined thread global, block sz: "<<idx<<","<<idy<<std::endl;
            idy = 0;
        }
    }
    for (unsigned idz =0; idz<idy; ++idz){
        thrs[idz].join();
        std::cout<< "Joined thread global, current, block sz: "<<idx<<","<<idz<<","<<idy<<std::endl;
    }
    //for (auto& th : thrs){
        //th.join();
    //} //finish up joining after jumping out

    //for (auto & g: gs) g.train();

    return 0;
};
