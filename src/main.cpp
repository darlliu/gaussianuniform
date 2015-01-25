#include "gucompute.h"

int main(int argc, char** argv){

    gurunner g;
    g.load("test_in.txt");
    g.train();
    return 0;
};
