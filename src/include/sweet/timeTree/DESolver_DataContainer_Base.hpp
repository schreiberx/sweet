/*
 * DESolver_DataContainer_Base.hpp
 */

#ifndef SRC_PROGRAMS_SIMDATA_DATACONTAINERBASE_HPP_
#define SRC_PROGRAMS_SIMDATA_DATACONTAINERBASE_HPP_


namespace sweet
{

class DESolver_DataContainer_Base
{
public:
	DESolver_DataContainer_Base()
	{
	}

public:
	virtual
	void swap(DESolver_DataContainer_Base &i_U) = 0;

public:
	virtual
	void op_setZero() = 0;

public:
	virtual
	void op_setVector(
			const DESolver_DataContainer_Base &i_a
		) = 0;

public:
	virtual
	void op_addVector(
			const DESolver_DataContainer_Base &i_a
		) = 0;

public:
	virtual
	void op_setVectorPlusVector(
			const DESolver_DataContainer_Base &i_a,
			const DESolver_DataContainer_Base &i_b
	) = 0;

public:
	virtual
	void op_setVectorPlusScalarMulVector(
			const DESolver_DataContainer_Base &i_a,
			double i_scalar,
			const DESolver_DataContainer_Base &i_b
		) = 0;

public:
	virtual
	void op_addScalarMulVector(
			double i_scalar,
			const DESolver_DataContainer_Base &i_a
		) = 0;
public:
	virtual
	void op_mulScalar(
			double i_scalar
		) = 0;

public:
	virtual
	void setup_like(
			const DESolver_DataContainer_Base &i_a
	) = 0;

	virtual
	DESolver_DataContainer_Base* getInstanceNew() const = 0;

	virtual
	void clear() = 0;

	virtual
	~DESolver_DataContainer_Base() {}
};

}

#endif
