

# Homebrew

Go to `https://brew.sh` and execute the installer script to install Homebrew.

Follow the instructions including how to add homebrew to your PATH.


# Homebrew packages


## Default packages

Install some required packages:
```
$ brew install wget
$ brew install cmake
```

## GCC

We also want to install a new gcc compiler:
```
$ brew install gcc@12
```

After this, the default compiler ```gcc``` will still point to the llvm compiler (not the GNU compiler!), but the gcc compiler is available as gcc-12.

## Alternative to miniconda

Next, we get a replacement for miniconda
```
$ brew install python@3.9
```

Start new terminal to load python environment and then install the following python packages:

```
$ pip3 install matplotlib numpy scipy
```

Again, start a new shell just to make sure that things are working right.

# Cloning SWEET
```
$ git clone https://github.com/schreiberx/sweet.git
```


# Setting up software for SWEET

## Change to SWEET directory
```
$ cd sweet
```

## Use bash! SWEET scripts are made to be used with 'bash'
```
$ bash
```

## Activate environment
```
$ source ./activate.sh
```


## Follow standard INSTALL instructions

E.g., install fftw, shtns, etc.

