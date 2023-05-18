/*
 * MyDataContainer.hpp
 */

#ifndef SRC_PROGRAMS_SIMDATA_MYDATACONTAINER_HPP_
#define SRC_PROGRAMS_SIMDATA_MYDATACONTAINER_HPP_

#include <algorithm>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/timeNew/PDESolver_DataContainer_Base.hpp>



class MyDataContainer :
	public sweet::PDESolver_DataContainer_Base
{
public:
	// Doesn't increase the size of the class
	static constexpr int Ndofs = 3;

	sweet::SphereData_Spectral data[Ndofs];

	// Note, that these references don't increase the size of the class
	sweet::SphereData_Spectral& phi = data[0];
	sweet::SphereData_Spectral& vrt = data[1];
	sweet::SphereData_Spectral& div = data[2];

public:
	MyDataContainer()
	{
	}

public:
	~MyDataContainer() override
	{
		clear();
	}

private:
	static inline
	MyDataContainer& cast(sweet::PDESolver_DataContainer_Base &i_U)
	{
		return static_cast<MyDataContainer&>(i_U);
	}

private:
	static inline
	const MyDataContainer& cast(const sweet::PDESolver_DataContainer_Base &i_U)
	{
		return static_cast<const MyDataContainer&>(i_U);
	}

public:
	void swap(
		PDESolver_DataContainer_Base &i_U
	) override
	{
		for (int i = 0; i < Ndofs; i++)
			data[i].swap(cast(i_U).data[i]);
	}

public:
	void setup_like(
		const sweet::PDESolver_DataContainer_Base &i_a
	) override
	{
		const MyDataContainer &i_d = static_cast<const MyDataContainer&>(i_a);

		for (int i = 0; i < Ndofs; i++)
			data[i].setup(i_d.data[i].sphereDataConfig);
	}

	PDESolver_DataContainer_Base* getNewInstance() const override
	{
		MyDataContainer *retval = new MyDataContainer;
		retval->setup_like(*this);

		return retval;
	}

public:
	void setup(
		const sweet::SphereData_Config *i_sphereData_Config
	)
	{
		for (int i = 0; i < Ndofs; i++)
			data[i].setup(i_sphereData_Config);
	}

public:
	void clear() override
	{
		for (int i = 0; i < Ndofs; i++)
			data[i].clear();
	}

public:
	void op_setVectorPlusVector(
			const sweet::PDESolver_DataContainer_Base &i_a,
			const sweet::PDESolver_DataContainer_Base &i_b
	) override
	{
		for (int i = 0; i < Ndofs; i++)
			data[i] = cast(i_a).data[i] + cast(i_b).data[i];
	}

public:
	void op_setVectorPlusScalarMulVector(
			const sweet::PDESolver_DataContainer_Base &i_a,
			double i_scalar,
			const sweet::PDESolver_DataContainer_Base &i_b
	) override
	{
		for (int i = 0; i < Ndofs; i++)
			data[i] = cast(i_a).data[i] + i_scalar*cast(i_b).data[i];
	}

public:
	void op_addScalarMulVector(
			double i_scalar,
			const sweet::PDESolver_DataContainer_Base &i_a
		) override
	{
		for (int i = 0; i < Ndofs; i++)
			data[i] += i_scalar*cast(i_a).data[i];
	}

};


#endif /* SRC_PROGRAMS_SIMDATA_MYDATACONTAINER_HPP_ */
