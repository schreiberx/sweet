
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
sweet::PlaneData_PhysicalComplex operator*(
		const double i_value,
		const sweet::PlaneData_PhysicalComplex &i_array_data
)
{
	return ((sweet::PlaneData_PhysicalComplex&)i_array_data)*i_value;
}


inline
static
sweet::PlaneData_PhysicalComplex operator*(
		const std::complex<double> &i_value,
		const sweet::PlaneData_PhysicalComplex &i_array_data
)
{
	return ((sweet::PlaneData_PhysicalComplex&)i_array_data)*i_value;
}
