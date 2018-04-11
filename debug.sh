rm test.so
python compile.py 1 2 3
g++ -g -std=c++11 test.cpp -o test.o
gdb test.o