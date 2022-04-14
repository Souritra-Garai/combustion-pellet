#ifndef __THERMAL_CONDUCTIVITY_PELLET__
#define __THERMAL_CONDUCTIVITY_PELLET__

namespace PelletThermalConductivity
{
	__device__ size_t max_iter = 1000;
	__device__ double tolerance = 0.001;


	__device__ double getThermalConductivityEMT(
		double v_2,
		double k_2,
		double k_1
	) {
		double v_1 = 1.0 - v_2;

		double k_1_k_2 = k_1 / k_2;

		double sqrt_D = sqrt(
			pow(k_1_k_2 * (3.0 * v_1 - 1.0), 2) +
			pow(2.0 - 3.0 * v_1, 2) +
			2.0 * k_1_k_2 * (2.0 + 9.0 * v_1 - 9.0 * v_1 * v_1)
		);

		return 0.25 * k_2 * (k_1_k_2 * (3.0 * v_1 - 1) + 2.0 - 3.0 * v_1 + sqrt_D);
	}

	__device__ double getThermalConductivityCC(
		double v_2,
		double k_2,
		double k_1
	) {
		double v_1 = 1.0 - v_2;

		double K_s = 1.0 / ((v_1 / k_1) + (v_2 / k_2));

		double K_p = k_1 * v_1 + k_2 * v_2;

		return (K_s / 2.0) * (sqrt(1.0 + 8.0 * K_p / K_s) - 1.0);
	}

	__device__ double getThermalConductivityME1(
		double v_2,
		double k_2,
		double k_1
	) {
		double v_1 = 1.0 - v_2;

		return (
			k_1 * v_1 * (2.0 * k_1 + k_2) + k_2 * v_2 * 3.0 * k_1
		) / (
			v_1 * (2.0 * k_1 + k_2) + v_2 * 3.0 * k_1
		);
	}


	__device__ double getThermalConductivityME2(
		double v_2,
		double k_2,
		double k_1
	) {
		double v_1 = 1.0 - v_2;

		return (
			k_2 * v_2 * (2.0 * k_2 + k_1) + k_1 * v_1 * 3.0 * k_2
		) / (
			v_2 * (2.0 * k_2 + k_1) + v_1 * 3.0 * k_2
		);
	}

	__device__ double getThermalConductivityME2(
		double v_1,
		double k_1,
		double v_2, 
		double k_2,
		double phi_11
	) {
		double f_1 = (1.0 - 2.0 * v_1 * phi_11) / 2.0;
		double f_2 = v_1 * phi_11 * 3.0 * k_2 / (2.0 * k_2 + k_1);

		return (f_1 * k_2 + f_2 * k_1) / (f_1 + f_2);
	}

	__device__ double getThermalConductivityEMT(
		double v_1,
		double k_1,
		double v_2, 
		double k_2,
		double phi_11
	) {
		double D = 
			(2.0 * k_1 - k_2) * v_1 * (1.0 - phi_11) + 
			(2.0 * k_2 - k_1) * (2.0 * v_2 + 2.0 * v_1 * phi_11 - 1.0) / 2.0
		;

		return 0.5 * (D + sqrt(D * D + 2.0 * k_1 * k_2));
	}

	__device__ double getThermalConductivityMEB(
		double v_2,
		double k_2,
		double k_1
	) {
		double v_1 = 1.0 - v_2;

		double phi_11_u = (2.0 * k_2 + k_1) / (2.0 * v_1 * (k_1 - k_2)) - 0.0000001;
		double phi_11_l = 0.0;
		
		double k_e_u = 
			getThermalConductivityEMT(v_1, k_1, v_2, k_2, phi_11_u) - 
			getThermalConductivityME2(v_1, k_1, v_2, k_2, phi_11_u);

		double k_e_l = 
			getThermalConductivityEMT(v_1, k_1, v_2, k_2, phi_11_l) - 
			getThermalConductivityME2(v_1, k_1, v_2, k_2, phi_11_l);

		unsigned int __iter = 0;

		double k_e_m, phi_11_m;

		while (abs(k_e_u - k_e_l) > tolerance && __iter < max_iter)
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

} // namespace PelletThermalConductivity

#endif