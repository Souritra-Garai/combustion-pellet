#include <thermo-physical-properties/Thermal_Conductivity_Pellet.hpp>

#include <math.h>

#define MAX_ITER 1000

template<typename real_t>
real_t getThermalConductivityCC(
	real_t v_2,
	real_t k_2,
	real_t k_1
) {
	real_t v_1 = 1.0 - v_2;

	real_t K_s = 1.0 / ((v_1 / k_1) + (v_2 / k_2));

	real_t K_p = k_1 * v_1 + k_2 * v_2;

	return (K_s / 2.0) * (sqrt(1.0 + 8.0 * K_p / K_s) - 1.0);
}

template<typename real_t>
real_t getThermalConductivityEMT(
	real_t v_2,
	real_t k_2,
	real_t k_1
) {
	real_t v_1 = 1.0 - v_2;

	real_t k_1_k_2 = k_1 / k_2;

	real_t sqrt_D = sqrt(
		pow(k_1_k_2 * (3.0 * v_1 - 1.0), 2) +
		pow(2.0 - 3.0 * v_1, 2) +
		2.0 * k_1_k_2 * (2.0 + 9.0 * v_1 - 9.0 * v_1 * v_1)
	);

	return 0.25 * k_2 * (k_1_k_2 * (3.0 * v_1 - 1) + 2.0 - 3.0 * v_1 + sqrt_D);
}

template<typename real_t>
real_t getThermalConductivityME1(
	real_t v_2,
	real_t k_2,
	real_t k_1
) {
	real_t v_1 = 1.0 - v_2;

	return (
		k_1 * v_1 * (2.0 * k_1 + k_2) + k_2 * v_2 * 3.0 * k_1
	) / (
		v_1 * (2.0 * k_1 + k_2) + v_2 * 3.0 * k_1
	);
}

template<typename real_t>
real_t getThermalConductivityME2(
	real_t v_2,
	real_t k_2,
	real_t k_1
) {
	real_t v_1 = 1.0 - v_2;

	return (
		k_2 * v_2 * (2.0 * k_2 + k_1) + k_1 * v_1 * 3.0 * k_2
	) / (
		v_2 * (2.0 * k_2 + k_1) + v_1 * 3.0 * k_2
	);
}

template<typename real_t>
real_t getThermalConductivityME2(
	real_t v_1,
	real_t k_1,
	real_t v_2, 
	real_t k_2,
	real_t phi_11
) {
	real_t f_1 = (1.0 - 2.0 * v_1 * phi_11) / 2.0;
	real_t f_2 = v_1 * phi_11 * 3.0 * k_2 / (2.0 * k_2 + k_1);

	return (f_1 * k_2 + f_2 * k_1) / (f_1 + f_2);
}

template<typename real_t>
real_t getThermalConductivityEMT(
	real_t v_1,
	real_t k_1,
	real_t v_2, 
	real_t k_2,
	real_t phi_11
) {
	real_t D = 
		(2.0 * k_1 - k_2) * v_1 * (1.0 - phi_11) + 
		(2.0 * k_2 - k_1) * (2.0 * v_2 + 2.0 * v_1 * phi_11 - 1.0) / 2.0
	;

	return 0.5 * (D + sqrt(D * D + 2.0 * k_1 * k_2));
}

template<typename real_t>
real_t getThermalConductivityMEB(
	real_t v_2,
	real_t k_2,
	real_t k_1
) {
	real_t v_1 = 1.0 - v_2;

	real_t phi_11_u = (2.0 * k_2 + k_1) / (2.0 * v_1 * (k_1 - k_2)) - 0.0000001;
	real_t phi_11_l = 0.0;
	
	real_t k_e_u = 
		getThermalConductivityEMT(v_1, k_1, v_2, k_2, phi_11_u) - 
		getThermalConductivityME2(v_1, k_1, v_2, k_2, phi_11_u);

	real_t k_e_l = 
		getThermalConductivityEMT(v_1, k_1, v_2, k_2, phi_11_l) - 
		getThermalConductivityME2(v_1, k_1, v_2, k_2, phi_11_l);

	unsigned int __iter = 0;

	real_t k_e_m, phi_11_m;

	while (abs(k_e_u - k_e_l) > 0.001 && __iter < MAX_ITER)
	{
		phi_11_m = 0.5 * (phi_11_u + phi_11_l);
		
		k_e_m = 
			getThermalConductivityEMT(v_1, k_1, v_2, k_2, phi_11_m) - 
			getThermalConductivityME2(v_1, k_1, v_2, k_2, phi_11_m);

		if (k_e_m == 0)
		
			break;

		else if (k_e_l * k_e_m < 0.0)
		{
			phi_11_u = phi_11_m;
			k_e_u = k_e_m;
		}

		else
		{
			phi_11_l = phi_11_m;
			k_e_l = k_e_m;
		}

		__iter++;
	}

	phi_11_m = 0.5 * (phi_11_u + phi_11_l);

	return	getThermalConductivityEMT(v_1, k_1, v_2, k_2, phi_11_m);
}

template long double	getThermalConductivityCC(	long double,	long double,	long double);
template double			getThermalConductivityCC(	double,			double,			double);
template float			getThermalConductivityCC(	float,			float,			float);

template long double	getThermalConductivityEMT(	long double,	long double,	long double);
template double			getThermalConductivityEMT(	double,			double,			double);
template float			getThermalConductivityEMT(	float,			float,			float);

template long double	getThermalConductivityME1(	long double,	long double,	long double);
template double			getThermalConductivityME1(	double,			double,			double);
template float			getThermalConductivityME1(	float,			float,			float);

template long double	getThermalConductivityME2(	long double,	long double,	long double);
template double			getThermalConductivityME2(	double,			double,			double);
template float			getThermalConductivityME2(	float,			float,			float);

template long double	getThermalConductivityMEB(	long double,	long double,	long double);
template double			getThermalConductivityMEB(	double,			double,			double);
template float			getThermalConductivityMEB(	float,			float,			float);