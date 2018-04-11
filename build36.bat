cls
del model_wrap.cpp _model.pyd model.pyc model.py
PATH = C:\Users\ROOZBEH\Downloads\swigwin-3.0.12;C:\Users\ROOZBEH\Downloads\mingw32\bin;C:\Python36
python compile.py
swig -c++ -python -o model_wrap.cpp model.i
g++ -shared -std=c++11 model_wrap.cpp -I C:\Python36\include -I C:\Python36\Lib\site-packages\numpy\core\include\ -LC:\Python36\libs\ -lpython36 -o _model.pyd
python plot_standalone.py
