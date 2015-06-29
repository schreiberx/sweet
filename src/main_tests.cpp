
#include "DataArray.hpp"
#include <unistd.h>


#define DIM 2


void runTests(
	std::size_t size_y,
	std::size_t size_x
)
{
	std::size_t size[2] = {size_x,size_y};

	DataArray<2> h(size);
	DataArray<2> hu(size);
	DataArray<2> hv(size);

	h.data_setall(10);
	hu.data_setall(0);
	hv.data_setall(0);

	int c = 1;
	for (std::size_t j = 0; j < size[1]; j++)
		for (std::size_t i = 0; i < size[0]; i++)
		{
			h.getDataRef(j, i) = c;
			c++;
		}

	// shift left test
	if (1)
	{
		DataArray<2> h_t(size);
		DataArray<2> u_t(size);
		DataArray<2> v_t(size);

		DataArray<2> op_shift_left(size);
		double shift_left_kernel[3][3] = {{0,0,0},{0,0,1},{0,0,0}};
		op_shift_left.setup_kernel(shift_left_kernel);

		std::cout << "H 0" << std::endl;
		std::cout << h << std::endl;

		h_t = op_shift_left(h);
		std::cout << "H 1" << std::endl;
		std::cout << h_t << std::endl;

		std::cout << h_t << std::endl;
		u_t = op_shift_left(h_t);
		std::cout << "H 2" << std::endl;
		std::cout << u_t << std::endl;

		h = op_shift_left(h);
		h = op_shift_left(h);
		std::cout << "H 2" << std::endl;
		std::cout << h << std::endl;

		return;
	}

	// shift up
	if (1)
	{
		DataArray<2> op_shift_left(size);
		double shift_left_kernel[3][3] = {{0,1,0},{0,0,0},{0,0,0}};
		op_shift_left.setup_kernel(shift_left_kernel);
		DataArray<2> shifted_result = op_shift_left(h);

		std::cout << "H" << std::endl;
		std::cout << h << std::endl;
		std::cout << "SHIFT UP" << std::endl;
		std::cout << shifted_result << std::endl;
	}

#if 0
	DataArray<2> test(size);
	double test_kernel[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
	//double test_kernel[5][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15},{16,17,18,19,20},{21,22,23,24,25}};
	test.setup_kernel(test_kernel);
	std::cout << test << std::endl;
	test = h;
	test.test_fftw();
	std::cout << test << std::endl;
	return;
#endif
}



int main(int i_argc, char *i_argv[])
{
	runTests(8, 8);
	return 1;
}
