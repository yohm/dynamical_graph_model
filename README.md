# dynamical_graph_model

Simulation code of the dynamical graph model proposed in the paper
Y. Murase et al., "A simple model for skewed species-lifetime distributions", New J. Phys. (2010).
http://iopscience.iop.org/article/10.1088/1367-2630/12/6/063021

# Usage

## Prerequisites

For analysis, ruby and gnuplot script are used.
Install these using your package-management system.
For Mac OS, run

```
brew install gnuplot
```

to install gnuplot. Ruby is installed by default.

## Build & Run

To build the simulation code, run

```
make
```

After the program was built, execute

```
./run.sh 0.2 1024 65536 1234
```

where the arguements specify "connection probability", "initial warm-up timesteps", "simulation timesteps", and the seed of random numbers, respectively.

You'll find the simulation data as well as png files of the plots.

# License

The MIT License (MIT)

Copyright (c) 2016 Yohsuke Murase

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

