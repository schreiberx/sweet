
inline
static
sweet::PlaneData_SpectralComplex operator*(
		const double i_value,
		const sweet::PlaneData_SpectralComplex &i_array_data
)
{
	return i_array_data*i_value;
}


inline
static
sweet::PlaneData_SpectralComplex operator*(
		const std::complex<double> &i_value,
		const sweet::PlaneData_SpectralComplex &i_array_data
)
{
	return i_array_data*i_value;
}


/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 *
 * Otherwise, we'd have to write it as arrayData+1.5
 *
 */

inline
static
sweet::PlaneData_SpectralComplex operator+(
		const double i_value,
		const sweet::PlaneData_SpectralComplex &i_array_data
)
{
	return ((sweet::PlaneData_SpectralComplex&)i_array_data)+i_value;
}

inline
static
sweet::PlaneData_SpectralComplex operator+(
		const std::complex<double> &i_value,
		const sweet::PlaneData_SpectralComplex &i_array_data
)
{
	return i_array_data+i_value;
}

inline
static
sweet::PlaneData_SpectralComplex operator-(
		const std::complex<double> &i_value,
		const sweet::PlaneData_SpectralComplex &i_array_data
)
{
	sweet::PlaneData_SpectralComplex out_plane_data(i_array_data.planeDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < i_array_data.planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_plane_data.spectral_space_data[idx] = -i_array_data.spectral_space_data[idx];

	out_plane_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

	return out_plane_data;

}
