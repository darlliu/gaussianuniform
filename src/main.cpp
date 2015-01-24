#include "gucompute.h"

int main(int argc, char** argv){

    gurunner g;
    g.load("test");
    g.run(1e-5, 50);
    return 0;
};
