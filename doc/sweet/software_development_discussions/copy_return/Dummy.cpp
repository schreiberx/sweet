#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <cstddef>
#include <Dummy.hpp>



Dummy2::Dummy2()
{
}


void Dummy::outputData()
{
	output("outputData()");
	for (int i = 0; i < data.size(); i++)
		std::cout << data[i] << std::endl;
}


Dummy::Dummy()
{
	output("Dummy()");
}

Dummy::Dummy(const Dummy &i_dummy)	:
	data(i_dummy.data)
{
	output("Dummy(&dummy)");

	std::stringstream buf;
	buf << "inh(" << i_dummy.id << ")";
	id = buf.str();
}

Dummy::Dummy(const Dummy2 &&i_dummy)
{
	output("Dummy(Dummy2 &&dummy)");
}

Dummy::Dummy(const Dummy &&i_dummy)	:
	data(i_dummy.data)
{
	output("Dummy(&&dummy)");

	std::stringstream buf;
	buf << "inh(" << i_dummy.id << ")";
	id = buf.str();
}

Dummy::Dummy(const std::string &i_id)	:
	id(i_id)
{
	output("Dummy(id)");

}

Dummy& Dummy::operator=(const Dummy &i_dummy)
{
	output("operator=(&dummy)");

	data = i_dummy.data;
	return *this;
}

Dummy& Dummy::operator=(const Dummy &&i_dummy)
{
	output("operator=(&&dummy)");

	data = i_dummy.data;
	return *this;
}

Dummy& Dummy::operator=(const Dummy2 &&i_dummy)
{
	output("operator=(Dummy2 &&dummy)");

	return *this;
}

Dummy Dummy::operator+(const Dummy &i_dummy)
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

Dummy::~Dummy()
{
	output("~Dummy()");
}

void Dummy::output(const std::string i_method_name)
{
	std::cout << "Class '" << id << "' method '" << i_method_name << "'" << std::endl;
}
