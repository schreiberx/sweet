#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <cstddef>


class Dummy2
{
public:
	Dummy2()
	{
	}
};


class Dummy
{
private:
	std::string id;
public:
	std::vector<int> data;

public:
	void outputData()
	{
		output("outputData()");
		for (int i = 0; i < data.size(); i++)
			std::cout << data[i] << std::endl;
	}


public:
	Dummy()
	{
		output("Dummy()");
	}

	Dummy(const Dummy &i_dummy)	:
		data(i_dummy.data)
	{
		output("Dummy(&dummy)");

		std::stringstream buf;
		buf << "inh(" << i_dummy.id << ")";
		id = buf.str();
	}

	Dummy(const Dummy2 &&i_dummy)
	{
		output("Dummy(Dummy2 &&dummy)");
	}

	__attribute__((noinline))
	Dummy(const Dummy &&i_dummy)	:
		data(i_dummy.data)
	{
		output("Dummy(&&dummy)");

		std::stringstream buf;
		buf << "inh(" << i_dummy.id << ")";
		id = buf.str();
	}

	Dummy(const std::string &i_id)	:
		id(i_id)
	{
		output("Dummy(id)");

	}

	Dummy& operator=(const Dummy &i_dummy)
	{
		output("operator=(&dummy)");

		data = i_dummy.data;
		return *this;
	}

	__attribute__((noinline))
	Dummy& operator=(const Dummy &&i_dummy)
	{
		output("operator=(&&dummy)");

		data = i_dummy.data;
		return *this;
	}

	__attribute__((noinline))
	Dummy& operator=(const Dummy2 &&i_dummy)
	{
		output("operator=(Dummy2 &&dummy)");

		return *this;
	}

	Dummy operator+(const Dummy &i_dummy)	const
	{
		Dummy ret;

		std::stringstream buf;
		buf << "(" << id << "+" << i_dummy.id << ")";
		ret.id = buf.str();

		ret.data.resize(data.size() + i_dummy.data.size());

		for (int i = 0; i < data.size(); i++)
			ret.data[i] = data[i];

		for (int i = 0; i < i_dummy.data.size(); i++)
			ret.data[i+data.size()] = i_dummy.data[i];

		return ret;
	}

	~Dummy()
	{
		output("~Dummy()");
	}

private:
	void output(const std::string i_method_name)
	{
		std::cout << "Class '" << id << "' method '" << i_method_name << "'" << std::endl;
	}
};


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
	std::cout << "* Dummy e = a+d;" << std::endl;
	Dummy e = a+d+a+d+b;
	e.outputData();

	std::cout << "**********************************************************" << std::endl;
	std::cout << "* Dummy f = Dummy2();" << std::endl;
	Dummy f = Dummy2();

	std::cout << "**********************************************************" << std::endl;
	std::cout << "* Dummy g;" << std::endl;
	std::cout << "* g = Dummy2();" << std::endl;
	Dummy g;
	g = Dummy2();

	std::cout << "**********************************************************" << std::endl;
	std::cout << "* return 0;" << std::endl;
	return 0;
}
