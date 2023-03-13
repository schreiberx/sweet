/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */


/**
 * operator to support operations such as:
 *
 * 1.5 * arrayData;
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
sweet::PlaneData_Physical operator*(
		const double i_value,
		const sweet::PlaneData_Physical &i_array_data
)
{
	return i_array_data*i_value;
}


/**
 * operator to support operations such as:
 *
 * 1.5 - arrayData;
 */
inline
static
sweet::PlaneData_Physical operator-(
		const double i_value,
		const sweet::PlaneData_Physical &i_array_data
)
{
	return i_array_data.operator_scalar_sub_this(i_value);
}


/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 */
inline
static
sweet::PlaneData_Physical operator+(
		const double i_value,
		const sweet::PlaneData_Physical &i_array_data
)
{
	return i_array_data+i_value;
}

