# CAC2021
Source code for paper presented in CAC2021.
## Install
This project uses `gnuplot`, `cmake`, `Eigen3` and `sciplot`. Go check them out if you don't have them locally installed.
### macOS
Use `brew` to install the missing packages.
```
$ brew install gnuplot cmake eigen
```
And git clone the [`sciplot`](https://sciplot.github.io) into the include path

## Usage
To compile and run the source code, please follow these steps
```
$ mkdir build && cd build 
$ cmake ..
$ make
$ ./cac
```

## [Optional]Data Visualization
Install [Veusz](https://veusz.github.io) and open [cac2021.vsz](cac2021.vsz) to plot the figure showed in the paper.
