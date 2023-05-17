/*
 * PDESolver_DataContainer_Base.hpp
 */

#ifndef SRC_PROGRAMS_SIMDATA_DATACONTAINERBASE_HPP_
#define SRC_PROGRAMS_SIMDATA_DATACONTAINERBASE_HPP_


namespace sweet
{

class PDESolver_DataContainer_Base
{
public:
	PDESolver_DataContainer_Base()
	{
	}

public:
	virtual
	void swap(PDESolver_DataContainer_Base &i_U) = 0;

public:
	virtual
	void op_setVectorPlusVector(
			const PDESolver_DataContainer_Base &i_a,
			const PDESolver_DataContainer_Base &i_b
	) = 0;

public:
	virtual
	void op_setVectorPlusScalarMulVector(
			const PDESolver_DataContainer_Base &i_a,
			double i_scalar,
			const PDESolver_DataContainer_Base &i_b
		) = 0;

public:
	virtual
	void op_addScalarMulVector(
			double i_scalar,
			const PDESolver_DataContainer_Base &i_a
		) = 0;

public:
	virtual
	void setup_like(
			const PDESolver_DataContainer_Base &i_a
	) = 0;

	virtual
	PDESolver_DataContainer_Base* getNewInstance() const = 0;

	virtual
	void clear() = 0;

	virtual
	~PDESolver_DataContainer_Base() {}
};

}

#endif
