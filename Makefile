all:
	@echo "This makefile is only a convenient handler for 'make clean', please use 'scons', see INSTALL"


clean:
	rm -rf build
	rm -rf /tmp/scons_build_*
	rm -rf python_mods/*.pyc
