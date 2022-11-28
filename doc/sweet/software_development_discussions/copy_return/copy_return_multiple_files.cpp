#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <cstddef>

#include <Dummy.hpp>

Dummy test()
{
	Dummy test("test");
	test.data.resize(2);
	test.data[0] = 3;
	test.data[1] = 4;

	return test;
}

int main()
{
	std::cout << "**********************************************************" << std::endl;
	std::cout << "* Dummy a(\"a\"); ..." << std::endl;
	Dummy a("a");
	a.data.resize(2);
	a.data[0] = 1;
	a.data[1] = 2;
	
	std::cout << "**********************************************************" << std::endl;
	std::cout << "* Dummy b(\"b\");" << std::endl;
	Dummy b("b");

	std::cout << "**********************************************************" << std::endl;
	std::cout << "* b = a;" << std::endl;
	b = a;

	std::cout << "**********************************************************" << std::endl;
	std::cout << "* WARNING: MOVE CONSTRUCTOR OPTIMIZED AWAY ***************" << std::endl;
	std::cout << "**********************************************************" << std::endl;
	std::cout << "* Dummy c = a;" << std::endl;
	Dummy c = a;

	std::cout << "**********************************************************" << std::endl;
	std::cout << "* WARNING: MOVE CONSTRUCTOR OPTIMIZED AWAY ***************" << std::endl;
	std::cout << "**********************************************************" << std::endl;
	std::cout << "* Dummy d = test();" << std::endl;
	Dummy d = test();
	d.outputData();

	std::cout << "**********************************************************" << std::endl;
	std::cout << "* WARNING: MOVE CONSTRUCTOR OPTIMIZED AWAY ***************" << std::endl;
	std::cout << "**********************************************************" << std::endl;
	std::cout << "* Dummy e = a+d+a+d+b;" << std::endl;
	Dummy e = a+d+a+d+b;
	e.outputData();

	std::cout << "**********************************************************" << std::endl;
	std::cout << "* Dummy f = Dummy2();" << std::endl;
	Dummy f = Dummy2();

	std::cout << "**********************************************************" << std::endl;
	std::cout << "* Dummy g;" << std::endl;
	std::cout << "* g = Dummy2();" << std::endl;
	Dummy g = Dummy2();

	std::cout << "**********************************************************" << std::endl;
	std::cout << "* return 0;" << std::endl;
	return 0;
}
