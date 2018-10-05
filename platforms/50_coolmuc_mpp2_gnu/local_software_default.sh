        
#PKGS+=("install_gcc7.1.sh")
PKGS+=("install_fftw3.sh")
#PKGS+=("install_eigen3.sh")
PKGS+=("install_lapack.sh")
#PKGS+=("install_openssl.sh")
PKGS+=("install_python3.sh")
PKGS+=("install_scons3.sh")
PKGS+=("install_shtns.sh")
#PKGS+=("install_shtns_python.sh")

# GIT support
# + First, install openssl
PKGS+=("install_openssl.sh")
# + Second, install curl which uses openssl
PKGS+=("install_curl.sh")
# + Third, install git which is based on curl
PKGS+=("install_git.sh")

