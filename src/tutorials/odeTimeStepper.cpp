/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/ErrorBase.hpp>

#include <sweet/timeNew/DESolver_Config_Base.hpp>

// Terms of differential equation
#include <sweet/timeNew/DESolver_DETerm_Base.hpp>
#include <sweet/timeNew/DESolver_DETerm_Registry.hpp>

#include <sweet/timeNew/DESolver_TimeStepping_Assemblation.hpp>
#include <sweet/timeNew/DESolver_TimeStepping_StringParser.hpp>
#include <sweet/timeNew/DESolver_TimeStepping_Tree.hpp>

// Terms of time stepper
#include <sweet/timeNew/DESolver_TimeStepper_Base.hpp>
#include <sweet/timeNew/DESolver_TimeStepper_Registry.hpp>
#include <sweet/timeNew/DESolver_TimeStepper_ExplicitRungeKutta.hpp>


#include <algorithm>


/**
 * We generate a simple polynomial of a particular degree
 * in order to test the right order of the time stepping methods.
 */
class SomePolynomial
{
	std::vector<double> alphas;

public:
	SomePolynomial()
	{
	}

	void setup(int i_degree)
	{
		alphas.resize(i_degree);

		for (std::size_t i = 0; i < alphas.size(); i++)
		{
			alphas[i] = 1.0/(double)(i+1.0);
			alphas[i] = i+1.0;
		}
	}

	/**
	 * Evaluate polynomial
	 */
	double eval(double i_x)
	{
		double acc = 0;
		double x_pow = 1.0;

		for (std::size_t i = 0; i < alphas.size(); i++)
		{
			acc += alphas[i]*x_pow;
			x_pow *= i_x;
		}

		return acc;
	}

	/**
	 * Evaluate derivate of polynomial
	 */
	double eval_deriv(double i_x)
	{
		double acc = 0;
		double x_pow = 1.0;

		for (std::size_t i = 1; i < alphas.size(); i++)
		{
			acc += alphas[i]*x_pow*(double)i;
			x_pow *= i_x;
		}

		return acc;
	}

	void print()
	{
		for (std::size_t i = 0; i < alphas.size(); i++)
		{
			std::cout << "alpha[" << i << "]: " << alphas[i] << std::endl;
		}
	}
};


/**
 * The main shack for this program to parse program arguments
 */
class ShackODETimeStepper :
		public sweet::ShackInterface
{
public:
	/// String of time stepping method
	std::string timeSteppingString;
	int polynomialDegree;

	ShackODETimeStepper()	:
		timeSteppingString(""),
		polynomialDegree(4)
	{

	}

	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "ShackOdeTimeStepper:" << std::endl;
		std::cout << "	--timestepping-method [string]	String of time stepping method" << std::endl;
		std::cout << "	--polynomial-degree [int]	Degree of polynomial to use for time stepping" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--timestepping-method", timeSteppingString);

		i_pa.getArgumentValueByKey("--polynomial-degree", polynomialDegree);

		if (i_pa.error.exists())
			return error.forwardWithPositiveReturn(i_pa.error);

		return true;

	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "ODE TIME STEPPER:" << std::endl;
		std::cout << " + timeSteppingString: " << timeSteppingString << std::endl;
		std::cout << " + polynomialDegree: " << polynomialDegree << std::endl;
		std::cout << std::endl;
	}
};



/*
 * Create a data container which allows storing just a vector of double values.
 *
 * This container will only be used with the interfaces declared in DESolver_DataContainer_Base
 */
class MyDataContainer :
	public sweet::DESolver_DataContainer_Base
{
public:
	std::vector<double> data;

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
	MyDataContainer& cast(sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<MyDataContainer&>(i_U);
	}

private:
	static inline
	const MyDataContainer& cast(const sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<const MyDataContainer&>(i_U);
	}

public:
	void swap(
		DESolver_DataContainer_Base &i_U
	) override
	{
		for (std::size_t i = 0; i < data.size(); i++)
			std::swap(data[i], cast(i_U).data[i]);
	}

	/*
	 * A simple setup which allocates as many Do
	 */
public:
	void setup(int i_size)
	{
		data.resize(i_size);
	}

public:
	void setup_like(
		const sweet::DESolver_DataContainer_Base &i_a
	) override
	{
		const MyDataContainer &i_d = static_cast<const MyDataContainer&>(i_a);

		data.resize(i_d.data.size());
	}

	DESolver_DataContainer_Base* getNewInstance() const override
	{
		MyDataContainer *retval = new MyDataContainer;
		retval->setup_like(*this);

		return retval;
	}

public:
	void clear() override
	{
		data.clear();
	}

public:
	void op_setVectorPlusVector(
			const sweet::DESolver_DataContainer_Base &i_a,
			const sweet::DESolver_DataContainer_Base &i_b
	) override
	{
		for (std::size_t i = 0; i < data.size(); i++)
			data[i] = cast(i_a).data[i] + cast(i_b).data[i];
	}

public:
	void op_setVectorPlusScalarMulVector(
			const sweet::DESolver_DataContainer_Base &i_a,
			double i_scalar,
			const sweet::DESolver_DataContainer_Base &i_b
	) override
	{
		for (std::size_t i = 0; i < data.size(); i++)
			data[i] = cast(i_a).data[i] + i_scalar*cast(i_b).data[i];
	}

public:
	void op_addScalarMulVector(
			double i_scalar,
			const sweet::DESolver_DataContainer_Base &i_a
		) override
	{
		for (std::size_t i = 0; i < data.size(); i++)
			data[i] += i_scalar*cast(i_a).data[i];
	}
};


/*
 * A special class which is forwarded to all
 *
 *  - time stepper instances and
 *  - DE term instances
 *
 * to set up data buffers and other things which are required.
 */
class MyDESolver_Config:
		public sweet::DESolver_Config_Base
{
public:
	/*
	 * Just a pointer to an existing data container
	 */
	MyDataContainer *myDataContainer;

	/*
	 * A polynomial for the DE Term "A"
	 */
	SomePolynomial *somePolynomial_A;

	/*
	 * A polynomial for the DE Term "B"
	 */
	SomePolynomial *somePolynomial_B;

	/*
	 * Return a new instance of a data container.
	 *
	 * This is what will be used by time steppers and/or DE term implementations
	 */
	sweet::DESolver_DataContainer_Base* getNewDataContainerInstance(
			int i_id = -1
	) const override
	{
		MyDataContainer *retval = new MyDataContainer;
		retval->setup_like(*myDataContainer);

		return retval;
	}
};


/**
 * This is one of the ODE terms we can solve for
 */
class MyODETerm_A	:
		public sweet::DESolver_DETerm_Base
{
public:
	sweet::ErrorBase error;

private:
	ShackODETimeStepper *shackODETimeStepper;

	SomePolynomial *_somePolynomial;
	double _dt;

public:
	MyODETerm_A()	:
		_dt(-1)
	{
	}

	~MyODETerm_A()	override
	{
	}


private:
	static inline
	MyDataContainer& cast(sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<MyDataContainer&>(i_U);
	}

	static inline
	const MyDataContainer& cast(const sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<const MyDataContainer&>(i_U);
	}

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override
	{
		shackODETimeStepper = io_shackDict->getAutoRegistration<ShackODETimeStepper>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

	const char* getImplementedPDETerm()	override
	{
		return "f";
	}

	std::shared_ptr<sweet::DESolver_DETerm_Base> getNewInstance() override
	{
		return std::shared_ptr<sweet::DESolver_DETerm_Base>(new MyODETerm_A);
	}

	const MyDESolver_Config& cast(const sweet::DESolver_Config_Base &i_config)
	{
		return static_cast<const MyDESolver_Config&>(i_config);
	}

	virtual
	bool setupDETermConfig(
		const sweet::DESolver_Config_Base &i_deTermConfig
	) override
	{
		const MyDESolver_Config& myConfig = cast(i_deTermConfig);

		_somePolynomial = myConfig.somePolynomial_A;
		assert(_somePolynomial != nullptr);
		_somePolynomial->setup(shackODETimeStepper->polynomialDegree);
		return true;
	}

	void setTimestepSize(double i_dt) override
	{
		_dt = i_dt;
	}

	/*
	 * Return the time tendencies of the PDE term
	 */
	void eval_tendencies(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_timeStamp
	)	override
	{
		//const MyDataContainer &i = cast(i_U);
		MyDataContainer &o = cast(o_U);

		o.data[0] = _somePolynomial->eval_deriv(i_timeStamp);
	}
};





class ODETimeStepper
{
public:
	sweet::ErrorBase error;

private:
	sweet::ShackProgArgDictionary shackProgArgDict;
	ShackODETimeStepper *shackODETimeStepper;

	typedef double T;

	MyDataContainer U;
	MyDataContainer U_tmp;

public:
	double timeStamp;
	double maxTimeStamp;
	double timeStepSize;

	SomePolynomial somePolynomial_A;
	SomePolynomial somePolynomial_B;

public:
	ODETimeStepper()	:
		shackODETimeStepper(nullptr)
	{
	}

	bool setup_1_shacks(int argc, char *argv[])
	{
		// setup shacks
		shackProgArgDict.setup(argc, argv);
		shackODETimeStepper = shackProgArgDict.getAutoRegistration<ShackODETimeStepper>();
		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	bool setup_2_config()
	{
		timeStamp = 0;
		timeStepSize = 0.25;

		return true;
	}

	bool setup_3_data()
	{
		U.data.resize(1);

		U_tmp.setup_like(U);


		return true;
	}

	std::shared_ptr<sweet::DESolver_TimeStepper_Base> timeStepper;

	bool setup_4_timestepper()
	{
		/*
		 * Setup time stepping string parser and parse it
		 */
		sweet::DESolver_TimeSteppingStringParser tsStringParser;
		sweet::DESolver_TimeStepping_Tree tsTree;

		tsStringParser.genTimeSteppingTree(
				shackODETimeStepper->timeSteppingString,
				tsTree
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(tsStringParser);
		tsTree.print();

		/*
		 * Register PDE Terms
		 */
		sweet::DESolver_DETerm_Registry pdeTerm_registry;
		pdeTerm_registry.registerPDETerm<MyODETerm_A>();

		/*
		 * Register time steppers
		 */
		sweet::DESolver_TimeStepper_Registry timeStepper_registry;
		timeStepper_registry.registerTimeStepper<sweet::DESolver_TimeStepper_ExplicitRungeKutta>();

		/*
		 * Ready to assemble time stepper
		 */
		sweet::DESolver_TimeStepping_Assemblation tssa;
		tssa.setup(pdeTerm_registry, timeStepper_registry);
		tssa.assembleTimeStepperByTree(
			tsTree,
			timeStepper
		);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(tssa);

		/*
		 * Call the shack Registration
		 */
		timeStepper->shackRegistration(&shackProgArgDict);

		/*
		 * final configuration of time steppers and DE Terms
		 */
		MyDESolver_Config deConfig;
		deConfig.myDataContainer = &U;
		deConfig.somePolynomial_A = &somePolynomial_A;
		deConfig.somePolynomial_B = &somePolynomial_B;

		timeStepper->setupConfig(deConfig);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timeStepper);

		/*
		 * Set time step size
		 */
		timeStepper->setTimeStepSize(timeStepSize);

		return true;
	}


	bool setup_5_initialConditions()
	{
		U.data[0] = somePolynomial_A.eval(0);

		somePolynomial_A.print();
		return true;
	}


	void clear()
	{
	}

	bool doTimestep()
	{
		timeStepper->eval_timeIntegration(U, U_tmp, timeStamp);
		U.swap(U_tmp);

		timeStamp += timeStepSize;

		if (std::abs(timeStamp - maxTimeStamp) < 1e-10)
			return false;

		return true;
	}

	double getValue()
	{
		return U.data[0];
	}

	double getError()
	{
		return std::abs(U.data[0] - somePolynomial_A.eval(timeStamp));
	}
};



int main(int argc, char *argv[])
{
	ODETimeStepper odeTimeStepper;

	odeTimeStepper.setup_1_shacks(argc, argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(odeTimeStepper);

	odeTimeStepper.setup_2_config();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(odeTimeStepper);

	odeTimeStepper.setup_3_data();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(odeTimeStepper);

	odeTimeStepper.setup_4_timestepper();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(odeTimeStepper);

	odeTimeStepper.setup_5_initialConditions();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(odeTimeStepper);

	odeTimeStepper.maxTimeStamp = 10.0;

	bool stop = false;
	while (true)
	{
		std::cout << odeTimeStepper.timeStamp << ", value=" << odeTimeStepper.getValue() << ", error=" << odeTimeStepper.getError() << std::endl;

		if (stop)
			break;

		if (!odeTimeStepper.doTimestep())
			stop = true;
	}

	return 0;
}
