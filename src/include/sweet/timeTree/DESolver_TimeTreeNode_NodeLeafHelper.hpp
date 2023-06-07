#ifndef SRC_SWEET_TIMETREE_TIMESTEPPER_NODE_LEAF_HELPER_HPP_
#define SRC_SWEET_TIMETREE_TIMESTEPPER_NODE_LEAF_HELPER_HPP_

#include "DESolver_TimeTreeNode_Base.hpp"


namespace sweet
{

/*!
 * Helper class for interior tree node
 *
 * It provides default member variables which are often used
 */
template <typename MainInteriorNodeClass>
class DESolver_TimeTreeNode_NodeLeafHelper	:
		public DESolver_TimeTreeNode_Base
{
protected:
	double _timestepSize;
	double &_dt = _timestepSize;

	// Number of stages to allocate buffers
	std::vector<sweet::DESolver_DataContainer_Base*> _tmpDataContainer;


public:
	DESolver_TimeTreeNode_NodeLeafHelper():
		_timestepSize(-1)
	{
	}


	/*!
	 * Copy constructor
	 *
	 * Simply copy the raw data over here
	 */
	DESolver_TimeTreeNode_NodeLeafHelper(
		const DESolver_TimeTreeNode_NodeLeafHelper &i_src	///<! Source node to copy from
	)	:
		DESolver_TimeTreeNode_Base(i_src)
	{
		_timestepSize = i_src._timestepSize;

		_tmpDataContainer.resize(_tmpDataContainer.size());
		for (std::size_t i = 0; i < i_src._tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_src._tmpDataContainer[i]->getNewDataContainer();
	}

	virtual
	~DESolver_TimeTreeNode_NodeLeafHelper()
	{
		clear();
	}


	inline
	void setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;
	}

	void clear() override
	{
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			delete _tmpDataContainer[i];

		_tmpDataContainer.clear();
	}

#if 0
	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceNew()	override
	{
#if 1
		MainInteriorNodeClass *m = new MainInteriorNodeClass(static_cast<MainInteriorNodeClass&>(*this));
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(m);
#else
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new MainInteriorNodeClass);
#endif
	}
#endif

	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceCopy()	override
	{
		MainInteriorNodeClass *m = new MainInteriorNodeClass(static_cast<MainInteriorNodeClass&>(*this));
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(m);
	}
};

}

#endif
