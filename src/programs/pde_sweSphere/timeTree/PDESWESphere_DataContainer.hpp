/*
 * PDESWESphere_DataContainer.hpp
 */

#ifndef SRC_PROGRAMS_SIMDATA_MYDATACONTAINER_HPP_
#define SRC_PROGRAMS_SIMDATA_MYDATACONTAINER_HPP_

#include <utility>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>



class PDESWESphere_DataContainer :
	public sweet::DESolver_DataContainer_Base
{
public:
	// How many number of DoF arrays
	static constexpr int Ndofs = 3;

	sweet::SphereData_Spectral data[Ndofs];

	// Note, that these references don't increase the size of the class
	sweet::SphereData_Spectral& phi_pert = data[0];
	sweet::SphereData_Spectral& div = data[2];
	sweet::SphereData_Spectral& vrt = data[1];

public:
	PDESWESphere_DataContainer()
	{
	}

public:
	~PDESWESphere_DataContainer() override
	{
		clear();
	}

private:
	static inline
	PDESWESphere_DataContainer& cast(sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<PDESWESphere_DataContainer&>(i_U);
	}

private:
	static inline
	const PDESWESphere_DataContainer& cast(const sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<const PDESWESphere_DataContainer&>(i_U);
	}

public:
	void swap(
		DESolver_DataContainer_Base &i_U
	) override
	{
		for (int i = 0; i < Ndofs; i++)
			data[i].swap(cast(i_U).data[i]);
	}

public:
	void setup_like(
		const sweet::DESolver_DataContainer_Base &i_a
	) override
	{
		const PDESWESphere_DataContainer &i_d = static_cast<const PDESWESphere_DataContainer&>(i_a);

		for (int i = 0; i < Ndofs; i++)
			data[i].setup(i_d.data[i].sphereDataConfig);
	}

	DESolver_DataContainer_Base* getInstanceNew() const override
	{
		PDESWESphere_DataContainer *retval = new PDESWESphere_DataContainer;
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
	void op_setVector(
			const sweet::DESolver_DataContainer_Base &i_a
	) override
	{
		const PDESWESphere_DataContainer &i_A = cast(i_a);
		for (int i = 0; i < Ndofs; i++)
			data[i] = i_A.data[i];
	}

public:
	void op_addVector(
			const sweet::DESolver_DataContainer_Base &i_a
	) override
	{
		const PDESWESphere_DataContainer &i_A = cast(i_a);
		for (int i = 0; i < Ndofs; i++)
			data[i] += i_A.data[i];
	}

public:
	void op_subVector(
			const sweet::DESolver_DataContainer_Base &i_a
	) override
	{
		const PDESWESphere_DataContainer &i_A = cast(i_a);
		for (int i = 0; i < Ndofs; i++)
			data[i] -= i_A.data[i];
	}

public:
	void op_setZero() override
	{
		for (int i = 0; i < Ndofs; i++)
			data[i].spectral_set_zero();
	}

public:
	void op_setVectorPlusVector(
			const sweet::DESolver_DataContainer_Base &i_a,
			const sweet::DESolver_DataContainer_Base &i_b
	) override
	{
		const PDESWESphere_DataContainer &i_A = cast(i_a);
		const PDESWESphere_DataContainer &i_B = cast(i_b);
		for (int i = 0; i < Ndofs; i++)
			data[i] = i_A.data[i] + i_B.data[i];
	}

public:
	void op_setVectorPlusScalarMulVector(
			const sweet::DESolver_DataContainer_Base &i_a,
			double i_scalar,
			const sweet::DESolver_DataContainer_Base &i_b
	) override
	{
		const PDESWESphere_DataContainer &i_A = cast(i_a);
		const PDESWESphere_DataContainer &i_B = cast(i_b);
		for (int i = 0; i < Ndofs; i++)
			data[i] = i_A.data[i] + i_scalar*i_B.data[i];
	}

public:
	void op_addScalarMulVector(
			double i_scalar,
			const sweet::DESolver_DataContainer_Base &i_a
		) override
	{
		const PDESWESphere_DataContainer &i_A = cast(i_a);
		for (int i = 0; i < Ndofs; i++)
			data[i] += i_scalar*i_A.data[i];
	}

public:
	void op_mulScalar(
			double i_scalar
		) override
	{
		for (int i = 0; i < Ndofs; i++)
			data[i] *= i_scalar;
	}
};


#endif /* SRC_PROGRAMS_SIMDATA_MYDATACONTAINER_HPP_ */
