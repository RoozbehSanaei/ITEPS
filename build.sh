clear
set -o xtrace
rm _model.so model.py model.pyc model_wrap.cpp _model.pyd
python compile.py 1 2 3
swig -c++ -python -o model_wrap.cpp model.i
g++ -std=c++11 -shared  -I /usr/include/python2.7/ -fPIC model_wrap.cpp -o _model.so 
python2 plot_standalone.py
