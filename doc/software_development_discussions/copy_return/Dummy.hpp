#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <cstddef>



class Dummy2
{
public:
	Dummy2();
};

class Dummy
{
private:
	std::string id;
public:
	std::vector<int> data;

public:
	void outputData();


public:
	Dummy();

	Dummy(const Dummy &i_dummy);

	Dummy(const Dummy2 &&i_dummy);

	Dummy(const Dummy &&i_dummy);

	Dummy(const std::string &i_id);

	Dummy& operator=(const Dummy &i_dummy);

	Dummy& operator=(const Dummy &&i_dummy);

	Dummy& operator=(const Dummy2 &&i_dummy);

	Dummy operator+(const Dummy &i_dummy);

	~Dummy();

private:
	void output(const std::string i_method_name);
};
