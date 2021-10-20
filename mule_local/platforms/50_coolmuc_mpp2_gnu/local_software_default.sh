        
# First, install OpenSSL
# This is also important for python
PKGS+=("install_openssl.sh")

PKGS+=("install_anaconda.sh")
PKGS+=("install_scons3.sh")

#PKGS+=("install_gcc7.1.sh")
PKGS+=("install_fftw3.sh")
PKGS+=("install_shtns.sh")
#PKGS+=("install_eigen3.sh")

PKGS+=("install_lapack.sh")
#PKGS+=("install_shtns_python.sh")

# GIT support
#  + OpenSSL was updated before
#  + Next, CURL is updated
PKGS+=("install_curl.sh")
#  + Leave existing 'git' untouched, makes problems on this server
#    The existing git works with the new curl libraries
