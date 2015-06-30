SCONS_OPTS:=-Q -j4


all:	debug


release:	gnu_release
gnu_release:
	scons $(SCONS_OPTS) --compiler=gnu --mode=release


debug:	gnu_debug
gnu_debug:
	scons $(SCONS_OPTS) --compiler=gnu --mode=debug


clean:
	rm -rf build
	rm -rf /tmp/scons_build_*
