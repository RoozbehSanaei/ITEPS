%module model

%{
  #define SWIG_FILE_WITH_INIT
  #include "model.h"
%}

/* Include the NumPy typemaps library */
%include "numpy.i"


%init %{
  import_array();
%}

%apply (double ARGOUT_ARRAY2[ANY][ANY]) { (double results[results_size][output_length]) };



%include "model.h"