all:	copy_return_single_file copy_return_multiple_files

copy_return_multiple_files:
	g++ -std=c++11 copy_return_multiple_files.cpp -c -o copy_return_multiple_files.o -O0 -I./ -fno-elide-constructors
	g++ -std=c++11 Dummy.cpp -c -o Dummy.o -O0 -I./ -fno-elide-constructors
	g++ -o copy_return_multiple_files copy_return_multiple_files.o Dummy.o
	

copy_return_single_file:
	g++ -std=c++11 copy_return_single_file.cpp -o copy_return_single_file -O0 -fno-elide-constructors

clean:
	rm -f copy_return_single_file

	rm -f copy_return_multiple_files.o
	rm -f Dummy.o
	rm -f copy_return_multiple_files
