        
# First, install OpenSSL
# This is also important for python
PKGS+=("install_openssl.sh")

if true; then
	PKGS+=("install_python3.sh")
	PKGS+=("install_scons3.sh")

	#PKGS+=("install_gcc7.1.sh")
	PKGS+=("install_fftw3.sh")
	PKGS+=("install_shtns.sh")
	#PKGS+=("install_eigen3.sh")

	PKGS+=("install_lapack.sh")
	#PKGS+=("install_shtns_python.sh")
fi

# GIT support
# + install the certificates
# + Here it is important to use the existing curl version to download the package
#PKGS+=("install_cacerts.sh")
# + install openssl
PKGS+=("install_openssl.sh")
# + install curl which uses openssl
PKGS+=("install_curl.sh")
# + install git which is based on curl
PKGS+=("install_git.sh")

