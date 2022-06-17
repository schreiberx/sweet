all:
	@echo "This makefile is only a convenient handler for 'make clean', please use 'scons', see INSTALL"


clean:
	#scons -c
	rm -rf build
	rm -rf /tmp/$(USER)/scons_build_*
	rm -rf python_mods/*.pyc

	rm -rf shtns
	rm -rf shtns_fftw
	rm -rf sweet_fftw
