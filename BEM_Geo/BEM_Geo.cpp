
#define _USE_MATH_DEFINES
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <cassert>

// global constants
double R = 6378000.0; // Earth radius [m]
double r_remaining = 500000; // [m]
double GM = 398600500000000.0; // [m^3 * s^(-1)]
double uExact = GM / R; // exact surface potential
double qExact = GM / (R * R); // surface acceleration

// const char* dataFilename = "BL-902.dat";
const char* dataFilename = "BL-3602.dat";
// const char* elemDataFilename = "elem_902.dat";
const char* elemDataFilename = "elem_3602.dat";
const int N = 3602; // dataset size;
const int Nt_max = 7;
const int NGauss = 7; // number of Gauss points per element

// Homogeneous coordinates of Gauss' points on a triangular element
const double etha_1 [NGauss] = {
	0.333333333, 0.79742699, 0.10128651, 0.10128651, 0.05971587, 0.47014206, 0.47014206
};
const double etha_2 [NGauss] = {
	0.333333333, 0.10128651, 0.79742699, 0.10128651, 0.47014206, 0.05971587, 0.47014206
};
const double etha_3 [NGauss] = {
	0.333333333, 0.10128651, 0.10128651, 0.79742699, 0.47014206, 0.47014206, 0.05971587
};

// Gauss' weights
const double w[NGauss] = {
	0.225000000, 0.12593918, 0.12593918, 0.12593918, 0.13239415, 0.13239415, 0.13239415
};

// ============================================================================================

// ============= Helper functions =============================================================

void printExactSol(std::string name, int printLim, bool inRow = true) {
	if (inRow) {
		std::cout << name << " = " << std::endl;
		for (int i = 0; i < printLim; i++) {
			std::cout << " " << uExact;
		}
		std::cout << "  ... ";
		for (int i = N - printLim; i < N; i++) {
			std::cout << " " << uExact;
		}
		std::cout << std::endl;
	}
	else {
		std::string offset = std::string((name + " = ").length() + 1, ' ');
		std::cout << name << " = " << std::endl;
		for (int i = 0; i < printLim; i++) {
			std::cout << offset << uExact << std::endl;
		}
		for (int i = 0; i < 3; i++) std::cout << offset << "  ." << std::endl;
		for (int i = N - printLim; i < N; i++) {
			std::cout << offset << uExact << std::endl;
		}
		std::cout << std::endl;
	}
}

void printArray1(std::string name, double* a, int printLim, bool inRow = true) {
	if (inRow) {
		std::cout << name << " = " << std::endl;
		for (int i = 0; i < printLim; i++) {
			std::cout << " " << a[i];
		}
		std::cout << "  ... ";
		for (int i = N - printLim; i < N; i++) {
			std::cout << " " << a[i];
		}
		std::cout << std::endl;
	}
	else {
		std::string offset = std::string((name + " = ").length() + 1, ' ');
		std::cout << name << " = " << std::endl;
		for (int i = 0; i < printLim; i++) {
			std::cout << offset << a[i] << std::endl;
		}
		for (int i = 0; i < 3; i++) std::cout << offset << "  ." << std::endl;
		for (int i = N - printLim; i < N; i++) {
			std::cout << offset << a[i] << std::endl;
		}
		std::cout << std::endl;
	}
}

void printArrayVector3(std::string name, double* vx, double* vy, double* vz, int printLim) {
	std::string offset = std::string((name + " = ").length() + 1, ' ');
	std::cout << name << " = " << std::endl;
	for (int i = 0; i < printLim; i++) {
		std::cout << offset << vx[i] << " " << vy[i] << " " << vz[i] << std::endl;
	}
	for (int i = 0; i < 3; i++) std::cout << offset << "  ." << std::endl;
	for (int i = N - printLim; i < N; i++) {
		std::cout << offset << vx[i] << " " << vy[i] << " " << vz[i] << std::endl;
	}
	std::cout << std::endl;
}

void printVertNeighborArray(std::string name, int* e, int printLim) {
	std::string offset = std::string((name + " = ").length() + 1, ' ');
	std::cout << name << " = " << std::endl;
	for (int i = 0; i < printLim; i++) {
		std::cout << offset << e[7 * i] << offset << " : " << offset;
		for (int j = 1; j < 7; j++) std::cout << e[7 * i + j] << " ";
		std::cout << std::endl;
	}
	for (int i = 0; i < 3; i++) std::cout << offset << "  ." << std::endl;
	for (int i = N - printLim; i < N; i++) {
		std::cout << offset << e[7 * i] << offset << " : " << offset;
		for (int j = 1; j < 7; j++) std::cout << e[7 * i + j] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void printArray2(std::string name, double** A, int printLim) {
	std::string offset = std::string((name + " = ").length() + 1, ' ');
	std::cout << name << " = " << std::endl;
	for (int i = 0; i < printLim; i++) {
		std::cout << offset;
		for (int j = 0; j < printLim; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << "  ...  ";
		for (int j = N - printLim; j < N; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
	for (int i = 0; i < 3; i++) std::cout << offset << "  ." << std::endl;
	for (int i = N - printLim; i < N; i++) {
		std::cout << offset;
		for (int j = 0; j < printLim; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << "  ...  ";
		for (int j = N - printLim; j < N; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void loadElemData(int* e) {
	printf("Loading elem data ... \n");

	std::fstream dataFile;
	dataFile.open(elemDataFilename, std::fstream::in);

	if (!dataFile.is_open()) {
		printf("Unable to open file %s\n", elemDataFilename);
	}
	else {
		printf("%s opened successfully\n", elemDataFilename);

		std::string line;
		int i = 0, j;

		while (std::getline(dataFile, line)) {
			std::vector<std::string> tokens;
			std::string s_delimiter = " ";
			size_t pos = 0;
			while (pos < 100) {
				pos = line.find(s_delimiter);
				tokens.push_back(line.substr(0, line.find(s_delimiter)));
				line = line.erase(0, pos + s_delimiter.length());
			}

			for (j = 0; j < tokens.size(); j++) e[7 * i + j] = (std::stoi(tokens[j]) - (j > 0 ? 1 : 0));

			i++;
		}

		dataFile.close();
	}
}

double vectorDot(double* a, double* b) {
	double result = 0.;
	for (int i = 0; i < N; i++) {
		result += a[i] * b[i];
	}
	return result;
}

double vectorNorm(double* a) {
	return sqrt(vectorDot(a, a));
}

// ==== Bi_CGSTAB solver ===============================================================

void Bi_CGSTAB_solve(double** A, double* b, double* x) {
	// ctrl. constants
	int maxIter = 100;
	double tol = 1e-6;

	// iter vectors
	double* x_curr = new double[N];
	double* x_next = new double[N];

	double* r_curr = new double[N];
	double* r_next = new double[N];

	double* rp0 = new double[N];

	double* p_curr = new double[N];
	double* p_next = new double[N];

	double* s = new double[N];

	double* tmp = new double[N];
	double* tmp1 = new double[N];

	// iter scalars
	double omega, alpha, beta, norm;

	// x0 = (1000,1000,...,1000)
#pragma omp parallel for
	for (int i = 0; i < N; i++) x_curr[i] = 1000.;
	// r0 = b - A x0
	// choose rp0 such that <r0, rp0> != 0
	// p0 = r0
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		r_curr[i] = b[i];
		for (int j = 0; j < N; j++) {
			r_curr[i] -= A[i][j] * x_curr[j];
		}
		rp0[i] = r_curr[i] + 100;
		p_curr[i] = r_curr[i];
	}
	std::cout << "==================================================" << std::endl;
	std::cout << "----------- Initializing Bi-CGSTAB Method --------" << std::endl;
	printArray2("systemMatrix", A, 4);
	printArray1("systemRhs", b, 5);
	printArray1("x0", x_curr, 2);
	printArray1("r0", r_curr, 5);
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "------------ Launching iterations ----------------" << std::endl;
	// begin iterations
	for (int k = 0; k < maxIter; k++) {
#pragma omp parallel
#pragma omp single
		std::cout << "::: iter : " << k << std::endl;

		// alpha[k] = <r[k], rp0> / <Ap[k], rp0>

		double num = 0.; double den = 0.;
#pragma omp parallel for reduction (+:num, den)
		for (int i = 0; i < N; i++) {
			tmp[i] = 0.;
			for (int j = 0; j < N; j++) {
				tmp[i] += A[i][j] * p_curr[j];
			}
			num += r_curr[i] * rp0[i];
			den += tmp[i] * rp0[i];
		}
		alpha = num / den;

		// s[k] = r[k] - alpha[k] * A p[k]
#pragma omp parallel for 
		for (int i = 0; i < N; i++) {
			s[i] = r_curr[i] - alpha * tmp[i];
		}

		norm = vectorNorm(s);
		std::cout << "||s|| = " << norm << std::endl;
		if (norm < tol) {
			// x[k + 1] = x[k] + alpha[k] * p[k]
#pragma omp parallel for
			for (int i = 0; i < N; i++) {
				x_next[i] = x_curr[i] + alpha * p_curr[i];
			}

			std::cout << "||s|| < tol = " << tol << ", exiting iterations" << std::endl;
			break;
		}

		// omega[k] = <A s[k], s[k]> / <A s[k], A s[k]>

		num = 0; den = 0;
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			tmp[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp[i] += A[i][j] * s[j];
			}
			num += tmp[i] * s[i];
			den += tmp[i] * tmp[i];
		}
		omega = num / den;

		// x[k + 1] = x[k] + alpha[k] * p[k] + omega[k] * s[k]
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			x_next[i] = x_curr[i] + alpha * p_curr[i] + omega * s[i];
		}

		// r[k + 1] = s[k] - omega[k] * A s[k]
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			r_next[i] = s[i] - omega * tmp[i];
		}

		norm = vectorNorm(r_next);
#pragma omp parallel
#pragma omp single
		std::cout << "||r[k + 1]|| = " << norm << std::endl;
		if (norm < tol) {
#pragma omp parallel
#pragma omp single
			std::cout << "||r[k + 1]|| < tol = " << tol << ", exiting iterations" << std::endl;
			break;
		}

		// beta[k] = (alpha[k] / omega[k]) * <r[k + 1], rp0> / <r[k], rp0>

		num = 0; den = 0;
#pragma omp parallel for reduction(+: num, den)
		for (int i = 0; i < N; i++) {
			num += r_next[i] * rp0[i];
			den += r_curr[i] * rp0[i];
		}

		beta = (alpha / omega) * num / den;

		// p[k + 1] = r[k + 1] + beta[k] * (p[k] - omega[k] * A p[k])
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			tmp[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp[i] += A[i][j] * p_curr[j];
			}
			p_next[i] = r_next[i] + beta * (p_curr[i] - omega * tmp[i]);
		}

		norm = fabs(vectorDot(r_next, rp0));
#pragma omp parallel
#pragma omp single
		std::cout << "|< r[k + 1], rp0 >| = " << norm << std::endl;
		if (norm < tol) {
			// rp0 = r[k + 1]; p[k + 1] = r[k + 1]
#pragma omp parallel for
			for (int i = 0; i < N; i++) {
				rp0[i] = r_next[i]; p_next[i] = r_next[i];
			}
		}
		// current = next
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			x_curr[i] = x_next[i];
			r_curr[i] = r_next[i];
			p_curr[i] = p_next[i];
		}
#pragma omp parallel
#pragma omp single
		std::cout << "===> finishing iter " << k << std::endl;
	}

	// result: x = x_next
#pragma omp parallel for
	for (int i = 0; i < N; i++) x[i] = x_next[i];

	// clean up
	delete[] x_curr; delete[] x_next;
	delete[] r_curr; delete[] r_next;
	delete[] p_curr; delete[] p_next;
	delete[] s; delete[] tmp; delete[] tmp1;
}


void writeData(double* B, double* L, double* u) {
	std::fstream dataOut;
	dataOut.open("data.dat", std::fstream::out);

	for (int i = 0; i < N; i++) {
		dataOut << B[i] << " " << L[i] << " " << u[i] << std::endl;
	}

	dataOut.close();
}

// ============================================================================================

// ============================================================================================

int main() {

	// ============= Load Files ======================================
	auto startLoad = std::chrono::high_resolution_clock::now();
	printf("Loading data ... \n");

	std::ifstream dataFile;
	dataFile.open(dataFilename);

	if (!dataFile.is_open()) {
		printf("Unable to open file %s\n", dataFilename);
		return 0;
	}

	printf("%s opened successfully\n", dataFilename);

	double* B = new double[N];
	double* L = new double[N];
	double* H = new double[N];

	double* q = new double[N];

	double dummy;

	for (int i = 0; i < N; i++) {
		dataFile >> B[i];
		dataFile >> L[i];
		dataFile >> H[i];
		dataFile >> q[i];
		dataFile >> dummy;
	}

	dataFile.close();

	printf("...point data loaded\n");

	// === Element data ===========

	dataFile.open(elemDataFilename);
	if (!dataFile.is_open()) {
		printf("Unable to open file %s\n", elemDataFilename);
		return 0;
	}

	printf("%s opened successfully\n", elemDataFilename);

	int* e = new int[(size_t)Nt_max * N];
	
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < Nt_max; j++) {
			dataFile >> e[Nt_max * i + j];
			if (j > 0) {
				e[Nt_max * i + j]--;
			}
		}
	}

	dataFile.close();

	printf("...elem data loaded\n");

	auto endLoad = std::chrono::high_resolution_clock::now();

	printf("Transforming point data to Cartesian coords:\n");

	double* x = new double[N];
	double* y = new double[N];
	double* z = new double[N];

	auto startTransform = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
	for (int i = 0; i < N; i++) {

		x[i] = (R + H[i]) * cos(B[i] * (M_PI / 180)) * cos(L[i] * (M_PI / 180));
		y[i] = (R + H[i]) * cos(B[i] * (M_PI / 180)) * sin(L[i] * (M_PI / 180));
		z[i] = (R + H[i]) * sin(B[i] * (M_PI / 180));

		q[i] = 0.00001 * q[i];
	}

	auto endTransform = std::chrono::high_resolution_clock::now();

	printArrayVector3("x", x, y, z, 5);
	printArray1("q", q, 5);
	printVertNeighborArray("e", e, 7);


	// ==============================================================================================================
	// ========== BOUNDARY ELEMENT METHOD ===========================================================================


	// Generating system matrices from point and elem data:

	auto startMatrixGen = std::chrono::high_resolution_clock::now();

	std::cout << "Pre-computing F and G helpers ... " << std::endl;

	// triangle Gauss pts
	double* tGauss_x = new double[(size_t)Nt_max * NGauss * N];
	double* tGauss_y = new double[(size_t)Nt_max * NGauss * N];
	double* tGauss_z = new double[(size_t)Nt_max * NGauss * N];

	// triangle edge helper vectors
	double* v_x = new double[(size_t)Nt_max * NGauss * N];
	double* v_y = new double[(size_t)Nt_max * NGauss * N];
	double* v_z = new double[(size_t)Nt_max * NGauss * N];

#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		int Ntri = e[Nt_max * i];

		for (int t = 0; t < Ntri; t++) {
			
			int tId = t + 1;
			int s = (t == Ntri - 1 ? 1 : tId + 1);

			int i1 = e[Nt_max * i + tId];
			int i2 = e[Nt_max * i + s];

			for (int k = 0; k < NGauss; k++) {
				int gauss_id = Nt_max * NGauss * i + NGauss * t + k;
				tGauss_x[gauss_id] = etha_1[k] * x[i] + etha_2[k] * x[i1] + etha_3[k] * x[i2];
				tGauss_y[gauss_id] = etha_1[k] * y[i] + etha_2[k] * y[i1] + etha_3[k] * y[i2];
				tGauss_z[gauss_id] = etha_1[k] * z[i] + etha_2[k] * z[i1] + etha_3[k] * z[i2];
			}

			v_x[Nt_max * i + t] = x[i] - x[e[Nt_max * i + tId]];
			v_y[Nt_max * i + t] = y[i] - y[e[Nt_max * i + tId]];
			v_z[Nt_max * i + t] = z[i] - z[e[Nt_max * i + tId]];
		}
	}

	double* nt_x = new double[(size_t)Nt_max * N];
	double* nt_y = new double[(size_t)Nt_max * N];
	double* nt_z = new double[(size_t)Nt_max * N];
	double* A_t = new double[(size_t)Nt_max * N];

#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		int Ntri = e[Nt_max * i];

		for (int t = 0; t < Ntri; t++) {
			int s = (t == Ntri - 1 ? 0 : t + 1);
			int Ti = Nt_max * i;

			nt_x[Ti + t] = v_y[Ti + t] * v_z[Ti + s] - v_z[Ti + t] * v_y[Ti + s];
			nt_y[Ti + t] = v_z[Ti + t] * v_x[Ti + s] - v_x[Ti + t] * v_z[Ti + s];
			nt_z[Ti + t] = v_x[Ti + t] * v_y[Ti + s] - v_y[Ti + t] * v_x[Ti + s];

			A_t[Ti + t] = 0.5 * sqrt(nt_x[Ti + t] * nt_x[Ti + t] + nt_y[Ti + t] * nt_y[Ti + t] + nt_z[Ti + t] * nt_z[Ti + t]);

			nt_x[Ti + t] /= (2 * A_t[Ti + t]);
			nt_y[Ti + t] /= (2 * A_t[Ti + t]);
			nt_z[Ti + t] /= (2 * A_t[Ti + t]);
		}
	}

	/*
	// Area test:
	double totArea = 0.0;
	for (int i = 0; i < N; i++) {
		int Ntri = e[Nt_max * i];
		for (int t = 0; t < Ntri; t++) {
			totArea += A_t[Nt_max * i + t];
		}
	}
	totArea /= 3;
	double sphereArea = 4 * M_PI * R * R;
	printf("\ntotArea  =  %.6lf \nsphereArea = %.6lf\ndiff = %.6lf %% \n\n", totArea, sphereArea, (sphereArea - totArea) / sphereArea * 100);
	*/

	std::cout << "... Pre-computing done!" << std::endl << std::endl;

	std::cout << "===> Filling F and rhs ..." << std::endl;

	// Filling F and rhs:

	double** F = new double* [N];
	double* rhs = new double [N];

#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		F[i] = new double[N];
		rhs[i] = 0.0;

		for (int j = 0; j < N; j++) {

			F[i][j] = 0.0;
			double G_ij = 0.0;

			if (i == j) { // SINGULAR elements
				double G_sum = 0.0;
				
				int Ti = Nt_max * i;
				int Ntri = e[Ti];

				for (int t = 0; t < Ntri; t++) {
					int s = (t == Ntri - 1 ? 0 : t + 1);
					

					double lx = x[e[Ti + t + 1]] - x[e[Ti + s + 1]];
					double ly = y[e[Ti + t + 1]] - y[e[Ti + s + 1]];
					double lz = z[e[Ti + t + 1]] - z[e[Ti + s + 1]];

					double lt_0 = sqrt(v_x[Ti + t] * v_x[Ti + t] + v_y[Ti + t] * v_y[Ti + t] + v_z[Ti + t] * v_z[Ti + t]);
					double lt_1 = sqrt(v_x[Ti + s] * v_x[Ti + s] + v_y[Ti + s] * v_y[Ti + s] + v_z[Ti + s] * v_z[Ti + s]);
					double lt_2 = sqrt(lx * lx + ly * ly + lz * lz);
					double alpha = acos((lt_1 * lt_1 + lt_0 * lt_0 - lt_2 * lt_2) / (2 * lt_1 * lt_0));
					double beta = acos((lt_0 * lt_0 + lt_2 * lt_2 - lt_1 * lt_1) / (2 * lt_0 * lt_2));

					G_sum += A_t[Ti + t] / lt_2 * log(tan(0.5 * (alpha + beta)) / tan(0.5 * beta));
				}

				G_ij = 1.0 / (4 * M_PI) * G_sum;
				rhs[i] += G_ij * q[i];
			}
			else { // REGULAR elements

				double rx = x[i] - x[j];
				double ry = y[i] - y[j];
				double rz = z[i] - z[j];
				double r_ij = sqrt(rx * rx + ry * ry + rz * rz);

				double F_sum = 0.0;
				double G_sum = 0.0;

				int Tj = Nt_max * j;
				int Ntri = e[Tj];

				for (int t = 0; t < Ntri; t++) {
					double K_ij = -(rx * nt_x[Tj + t] + ry * nt_y[Tj + t] + rz * nt_z[Tj + t]);

					double GaussSum_F = 0.0;
					double GaussSum_G = 0.0;

					for (int k = 0; k < NGauss; k++) {
						int gauss_id = Nt_max * NGauss * j + NGauss * t + k;
						double r_ijk_x = x[i] - tGauss_x[gauss_id];
						double r_ijk_y = y[i] - tGauss_y[gauss_id];
						double r_ijk_z = z[i] - tGauss_z[gauss_id];
						double r_ijk = sqrt(r_ijk_x * r_ijk_x + r_ijk_y * r_ijk_y + r_ijk_z * r_ijk_z);

						GaussSum_F += 1.0 / (r_ijk * r_ijk * r_ijk) * w[k] * etha_1[k];
						GaussSum_G += 1.0 / r_ijk * w[k] * etha_1[k];
					}

					F_sum += A_t[Tj + t] * K_ij * GaussSum_F;
					G_sum += A_t[Tj + t] * GaussSum_G;
				}

				F[i][j] = 1.0 / (4 * M_PI) * F_sum;

				G_ij = 1.0 / (4 * M_PI) * G_sum;
				rhs[i] += G_ij * q[j];
			}

		} // end j

		// Diagonal elements of F
		for (int j = 0; j < N; j++) {
			if (i != j) F[i][i] += F[i][j];
		}
		F[i][i] = 1.0 - F[i][i];
	}

	delete[] nt_x; delete[] nt_y; delete[] nt_z; delete[] A_t;
	delete[] tGauss_x; delete[] tGauss_y; delete[] tGauss_z;
	delete[] v_x; delete[] v_y; delete[] v_z;

	auto endMatrixGen = std::chrono::high_resolution_clock::now();


	printArray2("F", F, 4);
	printArray1("rhs", rhs, 4, false);

	// ============================================================================================================
	// ========== END BOUNDARY ELEM METHOD ========================================================================

	double* u = new double[N]; // potential solution

	// Bi-CGSTAB solve:
	auto startBi_CGSTAB = std::chrono::high_resolution_clock::now();

	Bi_CGSTAB_solve(F, rhs, u);

	auto endBi_CGSTAB = std::chrono::high_resolution_clock::now();
	// ========================= Solution done ===============================================

	// print solution
	printArray1("u", u, 4);

	// print exact solution
	printExactSol("uExact", 4);

	writeData(B, L, u);

	auto endWrite = std::chrono::high_resolution_clock::now();

	// ================================ Summary ============================================

	printf("\n============================================================================\n");
	printf("======================= Program summary =====================================\n");
	std::cout << "data size = " << N << std::endl;
	std::chrono::duration<double> elapsedTotal = (endWrite - startLoad);
	std::cout << "total runtime :    " << elapsedTotal.count() << " s" << std::endl;
	std::cout << "--------------------------------------------------------------------------" << std::endl;
	std::chrono::duration<double> elapsedLoad = (endLoad - startLoad);
	std::cout << "loading data file :    " << elapsedLoad.count() << " s" << std::endl;
	std::chrono::duration<double> elapsedTransform = (endTransform - startTransform);
	std::cout << "transforming data to Cartesian :      " << elapsedTransform.count() << " s" << std::endl;
	std::cout << "..........................................................................." << std::endl;
	std::chrono::duration<double> elapsedMatrixGen = (endMatrixGen - startMatrixGen);
	std::cout << "generating BEM system matrix :    " << elapsedMatrixGen.count() << " s" << std::endl;
	std::cout << ".... Bi-CGSTAB: .........................................................." << std::endl;
	std::chrono::duration<double> elapsedBi_CGSTAB = (endBi_CGSTAB - startBi_CGSTAB);
	std::cout << "Bi-CGSTAB solution :    " << elapsedBi_CGSTAB.count() << " s" << std::endl;
	std::cout << "--------------------------------------------------------------------------" << std::endl;


	// clean up
	delete[] x; delete[] y; delete[] z;
	delete[] q;

	for (int i = 0; i < N; i++) delete[] F[i];
	delete[] F;

	return 1;
}