#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_REXI_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_REXI_HPP_

#include <vector>
#include <string>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_NodeInteriorHelper.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/expIntegration/ExpFunction.hpp>
#include <sweet/expIntegration/REXI.hpp>
#include <sweet/expIntegration/REXICoefficients.hpp>
#include <sweet/expIntegration/ShackExpIntegration.hpp>


namespace sweet
{

class DESolver_TimeStepper_REXI	:
	public DESolver_TimeTreeNode_NodeInteriorHelper<DESolver_TimeStepper_REXI>
{
private:
	/*!
	 * Shack to exponential integration information including REXI information
	 */
	sweet::ShackExpIntegration *_shackExpIntegration;

	/*!
	 * String describing the exponential function ("phi0", "phi1" or ...)
	 */
	std::string _expFunctionString;

	/*!
	 * Should all REXI coefficients be
	 * - true: preallocated (typically fast) or
	 * - false: computed on-the-fly (typically slow)
	 */
	bool _rexiPreallocation;

	/*!
	 * REXI coefficients
	 */
	std::vector<std::complex<double>> _rexi_alphas;
	std::vector<std::complex<double>> _rexi_betas;
	std::complex<double> _rexi_gamma;

	/*!
	 * Number of REXI terms processed by this MPI rank
	 */
	std::size_t _num_local_rexi_terms;

	/*!
	 * Size of block to be processed locally by one particular thread (e.g. OMP thread)
	 */
	//std::size_t _num_local_rexi_terms_per_thread;

	/*!
	 * Number of local OpenMP threads
	 */
	int _num_local_rexi_parallel_threads;

	/*!
	 * Number of global OpenMP threads
	 */
	int _num_global_threads;


#if SWEET_MPI
	// MPI communicator
	MPI_Comm _mpi_comm;

	// number of mpi ranks to be used
	int _mpi_comm_rank;

	// MPI ranks
	int _mpi_comm_size;
#endif

public:
	DESolver_TimeStepper_REXI()	:
		_shackExpIntegration(nullptr),
		_rexiPreallocation(true),
		_num_local_rexi_terms(-1),
		//_num_local_rexi_terms_per_thread(-1),
		_num_local_rexi_parallel_threads(-1),
		_num_global_threads(-1)
	{
		setEvalAvailable(EVAL_EXPONENTIAL);
		setEvalAvailable(EVAL_INTEGRATION);
	}


	~DESolver_TimeStepper_REXI()
	{
		clear();
	}


	/*!
	 * This constructor is used for creating copies of this time stepper.
	 *
	 * This is required, e.g., for different phi functions where
	 * this class is duplicated for each phi function.
	 */
	DESolver_TimeStepper_REXI(
			const DESolver_TimeStepper_REXI &i_src
	)	:
		DESolver_TimeTreeNode_NodeInteriorHelper(i_src)
	{
		_rexiPreallocation = i_src._rexiPreallocation;
		_shackExpIntegration = i_src._shackExpIntegration;

		assert(_shackExpIntegration != nullptr);

		_num_local_rexi_terms = i_src._num_local_rexi_terms;
		//_num_local_rexi_terms_per_thread = i_src._num_local_rexi_terms_per_thread;
		_num_local_rexi_parallel_threads = i_src._num_local_rexi_parallel_threads;
		_num_global_threads = i_src._num_global_threads;
	}


	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("rexi");
		retval.push_back("REXI");
		return retval;
	}

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)	override
	{
		_shackExpIntegration = io_shackDict->getAutoRegistration<sweet::ShackExpIntegration>();

		return DESolver_TimeTreeNode_NodeInteriorHelper::shackRegistration(io_shackDict);
	}


	bool _setupArgumentInternals()
	{
		if (_timeTreeNodes.size() != 1)
			return error.set("One time node term needs to be given"+getNewLineDebugMessage());

		if (_expFunctionString == "")
			_expFunctionString = "phi0";

		return true;
	}


	virtual
	bool setupTreeNodeByFunction(
			std::shared_ptr<sweet::DESolver_TimeStepping_Tree::Function> &i_function,
			sweet::DESolver_TimeStepping_Assemblation &i_tsAssemblation
	)	override
	{
		for (auto iter = i_function->arguments.begin(); iter != i_function->arguments.end(); iter++)
		{
			sweet::DESolver_TimeStepping_Tree::Argument *a = iter->get();

			switch(a->argType)
			{
			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_FUNCTION:
			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_FUNCTION:
				if (_timeTreeNodes.size() == 1)
					return error.set("Only one term for exponential integration supported");

				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);

				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_VALUE:

				if (a->key == "preallocation")
				{
					a->getValue(_rexiPreallocation);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
					break;
				}
				return error.set("Key not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:
				if (_timeTreeNodes.size() != 0)
					return error.set("Only one term for exponential integration supported");

				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByName(
						a->value,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			default:
				SWEETError("Internal error");
				return error.set("Internal error");
			}
		}

		// Provide debug message in case that something goes wrong with the arguments
		setDebugMessage(i_function->getDebugMessage());
		return _setupArgumentInternals();
	}

	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	) override
	{
		if (i_key == "expIntegrationFunction")
		{
			_expFunctionString = i_value;
			return true;
		}

		return false;
	}

#if 0
	void _getWorkloadStartEnd(
			std::size_t &o_start,
			std::size_t &o_end,
			int i_local_thread_id
	)
	{
		std::size_t max_N = _rexi_alphas.size();

		#if SWEET_THREADING_TIME_REXI || SWEET_MPI

			#if SWEET_MPI
				int global_thread_id = i_local_thread_id + _num_local_rexi_parallel_threads*_mpi_comm_rank;
			#else
				int global_thread_id = i_local_thread_id;
			#endif

			assert(_num_local_rexi_terms >= 1);
			//assert(_num_local_rexi_terms_per_thread >= 0);
			assert(global_thread_id >= 0);

			o_start = std::min(max_N, _num_local_rexi_terms_per_thread*global_thread_id);
			o_end = std::min(max_N, o_start+_num_local_rexi_terms_per_thread);

		#else

			o_start = 0;
			o_end = max_N;

		#endif

#pragma omp critical
			std::cout << i_local_thread_id << ": " << o_start << " -> " << o_end << std::endl;

	}
#endif

	bool setupConfigAndGetTimeStepperEval(
		const sweet::DESolver_Config_Base &i_deTermConfig,
		EVAL_TYPES i_evalType,
		DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper
	) override
	{
		//rexiPreallocation
		assert(_timeTreeNodes.size() == 1);

#if SWEET_THREADING_TIME_REXI
		_num_local_rexi_parallel_threads = omp_get_max_threads();

		if (_num_local_rexi_parallel_threads == 0)
			SWEETError("FATAL ERROR: omp_get_max_threads == 0");
#else
		_num_local_rexi_parallel_threads = 1;
#endif

#if SWEET_MPI
		_mpi_comm = MPI_COMM_WORLD;	// TODO: Make me more flexible in future versions
		MPI_Comm_rank(_mpi_comm, &_mpi_comm_rank);
		MPI_Comm_size(_mpi_comm, &_mpi_comm_size);

		_num_global_threads = _num_local_rexi_parallel_threads * _mpi_comm_size;
#else
		_num_global_threads = _num_local_rexi_parallel_threads;
#endif

		/*
		 * Load REXI coefficients
		 */
		sweet::REXICoefficients<double> rexiCoefficients;
		sweet::REXI<> rexi;

		if (_shackExpIntegration->exp_method == "direct")
			return error.set("Direct exponential method is available with EXP() time tree function");

		rexi.load(
				_shackExpIntegration,
				_expFunctionString,
				rexiCoefficients,
				_shackExpIntegration->verbosity
		);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(rexi);

		_rexi_alphas.resize(rexiCoefficients.alphas.size());
		for (std::size_t n = 0; n < _rexi_alphas.size(); n++)
		{
			_rexi_alphas[n] = -rexiCoefficients.alphas[n];
			//_rexi_alphas[n] = rexiCoefficients.alphas[n];
		}

		_rexi_betas.resize(rexiCoefficients.betas.size());
		for (std::size_t n = 0; n < _rexi_betas.size(); n++)
		{
			_rexi_betas[n] = rexiCoefficients.betas[n];
		}

		_rexi_gamma = rexiCoefficients.gamma;

		assert(_rexi_alphas.size() > 0);

		/*
		 * Compute for which REXI coefficients we're responsible for
		 */
		std::size_t N = _rexi_alphas.size();
#if SWEET_MPI
		_num_local_rexi_terms = N/_mpi_comm_size;
		if (_num_local_rexi_terms*_mpi_comm_size != N)
			_num_local_rexi_terms++;
#else
		_num_local_rexi_terms = N;
#endif

#if 0
		_num_local_rexi_terms_per_thread = N/_num_global_threads;
		if (_num_local_rexi_terms_per_thread*_num_global_threads != N)
			_num_local_rexi_terms_per_thread++;
#endif
		/*
		 * Allocate data structures, time steppers, etc.
		 */
		_timeTreeNodes.resize(_num_local_rexi_terms);
		_evalFuns.resize(_num_local_rexi_terms);
		_tmpDataContainer.resize(_num_local_rexi_terms);

		/*
		 * Try to care about NUMA domains
		 */
#if SWEET_THREADING_TIME_REXI
		#pragma omp parallel for schedule(static,1) \
				default(none)	\
				shared(_timeTreeNodes,i_deTermConfig,_tmpDataContainer)
#endif
		for (std::size_t i = 0; i < _timeTreeNodes.size(); i++)
		{
			/*
			 * Initialize time tree node
			 */

			// skip 1st thread. We keep it in for equal split of for loop
			if (i != 0)
			{
				_timeTreeNodes[i] = _timeTreeNodes[0]->getInstanceCopy();
			}

			/*
			 * Setup evaluation
			 */
			_timeTreeNodes[i]->setupConfigAndGetTimeStepperEval(i_deTermConfig, EVAL_REXI_TERM, _evalFuns[i]);
#if SWEET_THREADING_TIME_REXI
			error.forward(_timeTreeNodes[i]->error);
#else
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
#endif

			assert(_evalFuns[i] != nullptr);


			/*
			 * REXI alpha and beta coefficients
			 */
			_timeTreeNodes[i]->setupByKeyValue(
					"rexiTermAlpha",
					_rexi_alphas[i]
				);

			_timeTreeNodes[i]->setupByKeyValue(
					"rexiTermBeta",
					_rexi_betas[i]
				);


			/*
			 * setup data container for output, use argument "1" to request complex data
			 */
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

			if (_rexiPreallocation)
				_timeTreeNodes[i]->setupByKeyValue("rexiTermPreallocation", "true");
			else
				_timeTreeNodes[i]->setupByKeyValue("rexiTermPreallocation", "false");

#if SWEET_THREADING_TIME_REXI
			error.forward(_timeTreeNodes[i]->error);
#else
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
#endif
		}

		if (error.exists())
			return false;

		if (_num_local_rexi_parallel_threads == 0)
			SWEETError("FATAL ERROR C: omp_get_max_threads == 0");

		// default setup
		DESolver_TimeTreeNode_Base::_helperGetTimeStepperEval(
				i_evalType,
				o_timeStepper
			);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

		return true;
	}


#if 0
	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceNew()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_REXI);
	}
#endif

	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_REXI(*this));
	}

	inline
	void setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		assert(_timeTreeNodes.size() == _num_local_rexi_terms);

#if SWEET_THREADING_TIME_REXI
		#pragma omp parallel for schedule(static,1) \
			default(none)		\
			shared(_timeTreeNodes, i_dt)
#endif
		for (std::size_t local_thread_id = 0; local_thread_id < _timeTreeNodes.size(); local_thread_id++)
		{
			_timeTreeNodes[local_thread_id]->setTimeStepSize(i_dt);
		}
	}




private:
	bool _eval_integration(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)	override
	{
		return _eval_exponential(i_U, o_U, i_simulationTime);
	}


	/*!
	 * We also provide an exponential time integration for this one
	 * in order to transparently support 'exponential' time integration for
	 * either DE terms themselves, EXP and also REXI evaluations.
	 */
private:
	bool _eval_exponential(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)	override
	{
		const sweet::DESolver_DataContainer_Base *i_U_real;

#if SWEET_MPI
		o_U.op_setVector(i_U);
		o_U.mpiBcast(_mpi_comm);
		i_U_real = &o_U;
#else
		i_U_real = &i_U;
#endif


#if SWEET_THREADING_TIME_REXI
		#pragma omp parallel for schedule(static,1) \
		default(none)								\
		shared(_timeTreeNodes,i_U_real,_tmpDataContainer,i_simulationTime)
#endif
		for (std::size_t i = 0; i < _timeTreeNodes.size(); i++)
		{
			evalTimeStepper(
					i,
					*i_U_real,
					*_tmpDataContainer[i],
					i_simulationTime
				);
		}

		/*
		 * REDUCE operation
		 */
#if SWEET_MPI
		/*
		 * Step 1) Reduce to first tmpDataContainer.
		 * Step 2) Call MPIReduce
		 */
		for (std::size_t i = 1; i < _timeTreeNodes.size(); i++)
		{
			_tmpDataContainer[0]->op_addVector(*_tmpDataContainer[i]);
		}

		o_U.mpiReduce(
				*_tmpDataContainer[0],
				_mpi_comm
			);

#else
		o_U.op_setZero();
		for (std::size_t i = 0; i < _timeTreeNodes.size(); i++)
		{
			o_U.op_addVector(*_tmpDataContainer[i]);
		}
#endif

		if (_rexi_gamma.real() != 0 || _rexi_gamma.imag() != 0)
		{
			SWEETError("TODO");
//			o_U.op_addScalarMulVector(i_scalar, i_U);
		}

		return true;
	}


	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "REXI(" << std::endl;
		std::cout << newPrefix << "  expFunctionString: '" << _expFunctionString << "'" << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}

#endif
