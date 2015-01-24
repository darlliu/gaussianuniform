#!/bin/bash
cd src
clang++  -I/home/yul13/bin/include -L/home/yul13/bin/lib gucompute.h gucompute.cpp main.cpp  -lz -lstdc++ -lrt -lpthread -lm -lc -O3 -std=c++1y
mv a.out ../bin/gucompute.exe
cd ..
