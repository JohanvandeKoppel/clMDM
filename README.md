# clMDM
This repository contains the OpenCL code that runs the Mussel disturbance model describing spatial self-organisation of mussel beds on rocky shores. The model has been published in:

`Guichard, F., P. M. Halpin, G. W. Allison, J. Lubchenco and B. A. Menge (2003). "Mussel disturbance dynamics: signatures of oceanographic forcing from local interactions." American Naturalist 161(6): 889-904.`

The model is a cellular automaton, which requires the uses a random number generator. I use a Mersenne Twister PRNG, which I found in [this](https://github.com/jj9854/oclFeldmanCousins) repository.

An video of the output of the model can be found on: https://www.youtube.com/watch?v=RzrqpuVluGg. This video has been produced by the code in this repository.

The use of graphics processors for simulating ecological models (using CUDA) at large spatial scales is described in:

`Van de Koppel, J., R. Gupta and C. Vuik (2011). "Scaling-up spatially-explicit ecological models using graphics processors." Ecological Modelling 222(17): 3011-3019.`

To run this model, read the explanation on the main page https://johanvandekoppel.github.io.
