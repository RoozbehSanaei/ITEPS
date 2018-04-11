cls
del model_wrap.cpp _model.pyd model.pyc model.py
PATH = C:\Users\ROOZBEH\Downloads\swigwin-3.0.12;C:\Users\ROOZBEH\Downloads\mingw32\bin;C:\Python27;C:\Python27\Scripts
python compile.py %1 %2 %3
swig -c++ -python -o model_wrap.cpp model.i
g++ -shared -std=c++11 model_wrap.cpp -I C:\Python27\include -I C:\Python36\Lib\site-packages\numpy\core\include\ -LC:\Python27\libs\ -lpython27 -o _model.pyd
REM python plot_standalone.py
