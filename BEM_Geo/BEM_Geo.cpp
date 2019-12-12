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
#include <mpi.h>

// global constants
double R = 6378000.0; // Earth radius [m]
double r_remaining = 500000; // [m]
double GM = 398600500000000.0; // [m^3 * s^(-1)]
double uExact = GM / R; // exact surface potential
double qExact = GM / (R * R); // surface acceleration

const char* dataFilename = "BL-902.dat";
// const char* dataFilename = "BL-3602.dat";
const char* elemDataFilename = "elem_902.dat";
// const char* elemDataFilename = "elem_3602.dat";
const int N = 902; // dataset size;
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

void printArray1FromTo(std::string name, double* a, int from, int to, int printLim, bool inRow = true) {
	if (inRow) {
		std::cout << name << "[" << from << ":" << to << "] = ";
		for (int i = from; i < from + printLim; i++) {
			std::cout << " " << a[i];
		}
		std::cout << "  ... ";
		for (int i = to - printLim; i < to; i++) {
			std::cout << " " << a[i];
		}
		std::cout << std::endl;
	}
	else {
		std::string offset = std::string((name + " = ").length() + 1, ' ');
		std::cout << name << "[" << from << ":" << to << "] = " << std::endl;
		for (int i = from; i < from + printLim; i++) {
			std::cout << offset << a[i] << std::endl;
		}
		for (int i = 0; i < 3; i++) std::cout << offset << "  ." << std::endl;
		for (int i = to - printLim; i < to; i++) {
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

int main(int argc, char** argv) {
	// point data

	double* B = new double[N + 10];
	double* L = new double[N + 10];
	double* H = new double[N + 10];

	double* x = new double[N + 10];
	double* y = new double[N + 10];
	double* z = new double[N + 10];

	double* q = new double[N + 10];
	double* d2U = new double[N + 10];

	// elem data

	int* e = new int[(size_t)Nt_max * N + 10];

	// system rhs

	double* rhs = new double[N + 10];

	// Bi-CGSTAB
	// ctrl. constants
	const int maxIter = 100;
	const double tol = 2e-5;

	// iter vectors
	double* x_curr = new double[N + 10];
	double* x_next = new double[N + 10];

	double* r_curr = new double[N + 10];
	double* r_next = new double[N + 10];

	double* rp0 = new double[N + 10];

	double* p_curr = new double[N + 10];
	double* p_next = new double[N + 10];

	double* s = new double[N + 10];

	double* tmp = new double[N + 10];

	// iter scalars
	double alpha, beta, omega;

	// MPI Vars:
	int nprocs, myrank;
	int istart, iend, nlocal = 0, nlast = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	// ============= File Load =========================================

	auto startLoad = std::chrono::high_resolution_clock::now();
	if (myrank == 0) {
		std::fstream dataFile;

		printf("Loading data ... \n");

		dataFile.open(dataFilename, std::fstream::in);

		if (!dataFile.is_open()) {
			printf("Unable to open file %s\n", dataFilename);
			return 0;
		}
		else {
			printf("%s opened successfully\n", dataFilename);

			std::string line;
			int i = 0;

			while (std::getline(dataFile, line)) {
				std::vector<std::string> tokens;
				std::string s_delimiter = " ";
				size_t pos = 0;
				while (pos < 100) {
					pos = line.find(s_delimiter);
					tokens.push_back(line.substr(0, line.find(s_delimiter)));
					line = line.erase(0, pos + s_delimiter.length());
				}

				B[i] = std::stod(tokens[0]);
				L[i] = std::stod(tokens[1]);
				H[i] = std::stod(tokens[2]);
				double Q = std::stod(tokens[3]);
				double D2U = std::stod(tokens[4]);

				x[i] = (R + H[i]) * cos(B[i] * M_PI / 180) * cos(L[i] * M_PI / 180);
				y[i] = (R + H[i]) * cos(B[i] * M_PI / 180) * sin(L[i] * M_PI / 180);
				z[i] = (R + H[i]) * sin(B[i] * M_PI / 180);

				q[i] = 0.00001 * Q;
				d2U[i] = D2U;

				i++;
			}

			dataFile.close();
		}

		printf("Loading elem data ... \n");

		dataFile.open(elemDataFilename, std::fstream::in);

		if (!dataFile.is_open()) {
			printf("Unable to open file %s\n", elemDataFilename);
			return 0;
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
	auto endLoad = std::chrono::high_resolution_clock::now();

	auto startPrint3 = std::chrono::high_resolution_clock::now();
	if (myrank == 0) printArrayVector3("x", x, y, z, 5);
	auto endPrint3 = std::chrono::high_resolution_clock::now();

	auto startPrint1 = std::chrono::high_resolution_clock::now();
	if (myrank == 0) printArray1("q", q, 5);
	auto endPrint1 = std::chrono::high_resolution_clock::now();
	if (myrank == 0) printArray1("d2U/dn2", d2U, 5);

	if (myrank == 0) printVertNeighborArray("e", e, 7);

	// ========== Files loaded ===========================================

	// array broadcasts
	MPI_Bcast(q, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(z, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(e, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// local/last packet sizes
	nlocal = (N / nprocs) + 1;
	nlast = N - (nprocs - 1) * nlocal;

	// packet indexing
	istart = nlocal * myrank;
	if (myrank == nprocs - 1)
		iend = N - 1;
	else
		iend = istart + nlocal - 1;

	std::cout << "p" << myrank << " range: " << istart << " - " << iend << ", nlocal: " << nlocal << ", nlast = " << nlast << std::endl;


	// ==============================================================================================================
	// ========== BOUNDARY ELEMENT METHOD ===========================================================================


	// Generating system matrices from point and elem data:

	double* nt_x = new double[(size_t)Nt_max * (N + 10)];
	double* nt_y = new double[(size_t)Nt_max * (N + 10)];
	double* nt_z = new double[(size_t)Nt_max * (N + 10)];
	double* A_t = new double[(size_t)Nt_max * (N + 10)];

	double* G_diag = new double[N + 10]; // diagonal (singular) elems of G-matrix

	double* r_ijk_x = new double[(size_t)Nt_max * NGauss * (N + 10)];
	double* r_ijk_y = new double[(size_t)Nt_max * NGauss * (N + 10)];
	double* r_ijk_z = new double[(size_t)Nt_max * NGauss * (N + 10)];

	auto startMatrixGen = std::chrono::high_resolution_clock::now();
	if (myrank == 0) {
		std::cout << "Pre-computing F and G helpers ... " << std::endl;

		int Ntri, j, t, k;
		int j1, j2;
		double t0x, t0y, t0z, t1x, t1y, t1z, t2x, t2y, t2z;
		double nx, ny, nz;
		double norm;

		double e0x, e0y, e0z, e1x, e1y, e1z, e2x, e2y, e2z; // triangle edges
		double l0_t, l1_t, l2_t, At;

		for (j = 0; j < N; j++) {

			Ntri = e[Nt_max * j];
			G_diag[j] = 0.0;

			for (t = 1; t <= Ntri; t++) {
				// adjacent vertex indices
				j1 = e[Nt_max * j + t];
				j2 = e[Nt_max * j + t % Ntri + 1];

				// ===================== Precompute normals ===========================
				// tri verts
				t0x = x[j];		t0y = y[j];		t0z = z[j];
				t1x = x[j1];	t1y = y[j1];	t1z = z[j1];
				t2x = x[j2];	t2y = y[j2];	t2z = z[j2];

				// tri normal
				nx = (t1y - t0y) * (t2z - t0z) - (t1z - t0z) * (t2y - t0y);
				ny = (t1z - t0z) * (t2x - t0x) - (t1x - t0x) * (t2z - t0z);
				nz = (t1x - t0x) * (t2y - t0y) - (t1y - t0y) * (t2x - t0x);

				norm = sqrt(nx * nx + ny * ny + nz * nz);
				assert(norm != 0.0);

				nx /= norm;
				ny /= norm;
				nz /= norm;

				nt_x[Nt_max * j + t] = nx;
				nt_y[Nt_max * j + t] = ny;
				nt_z[Nt_max * j + t] = nz;

				// ==============  Pre-compute G-diagonal elems ===============================
				// e0 = t0 -> t1
				e0x = x[j1] - x[j];
				e0y = y[j1] - y[j];
				e0z = z[j1] - z[j];
				// e1 = t1 -> t2
				e1x = x[j2] - x[j1];
				e1y = y[j2] - y[j1];
				e1z = z[j2] - z[j1];
				// e2 = t2 -> t0
				e2x = x[j] - x[j2];
				e2y = y[j] - y[j2];
				e2z = z[j] - z[j2];

				// edge lengths
				l0_t = sqrt(e0x * e0x + e0y * e0y + e0z * e0z);
				l1_t = sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
				l2_t = sqrt(e2x * e2x + e2y * e2y + e2z * e2z);

				assert(l0_t != 0.0 && l1_t != 0.0 && l2_t != 0.0);

				alpha = acos(-(e0x * e2x + e0y * e2y + e0z * e2z) / (l0_t * l2_t));
				beta = acos(-(e0x * e1x + e0y * e1y + e0z * e1z) / (l0_t * l1_t));

				assert(alpha != 0.0 && beta != 0.0);

				At = 0.5 * norm; // triangle area

				A_t[Nt_max * j + t] = At;
				G_diag[j] += At / l1_t * log(tan(0.5 * (alpha + beta)) / tan(0.5 * beta));

				// ====================== Pre-compute triangle Gauss' pts coords ========================
				for (k = 0; k < NGauss; k++) {
					double x_ijk = etha_1[k] * x[j] + etha_2[k] * x[j1] + etha_3[k] * x[j2];
					double y_ijk = etha_1[k] * y[j] + etha_2[k] * y[j1] + etha_3[k] * y[j2];
					double z_ijk = etha_1[k] * z[j] + etha_2[k] * z[j1] + etha_3[k] * z[j2];

					r_ijk_x[Nt_max * NGauss * j + NGauss * t + k] = x_ijk;
					r_ijk_y[Nt_max * NGauss * j + NGauss * t + k] = y_ijk;
					r_ijk_z[Nt_max * NGauss * j + NGauss * t + k] = z_ijk;
				}
			}
		}

		std::cout << "... Pre-computing done!" << std::endl << std::endl;
	}

	// array broadcasts
	MPI_Bcast(nt_x, Nt_max * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(nt_y, Nt_max * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(nt_z, Nt_max * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(A_t, Nt_max * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(G_diag, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(r_ijk_x, Nt_max * NGauss * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(r_ijk_y, Nt_max * NGauss * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(r_ijk_z, Nt_max * NGauss * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (myrank == 0) std::cout << "===> Filling F and rhs ..." << std::endl;

	// Filling F and rhs:

	double** F_local = new double* [nlocal];
	double* rhs_local = new double [nlocal];
	double* F_diag = new double[N];

	for (int i = 0; i < nlocal; i++) F_local[i] = new double[N];

	int iGlobal = 0;
	for (int i = 0; i < nlocal; i++) {
		iGlobal = i + istart;

		double rhs_sum = 0.0;

		for (int j = 0; j < N; j++) {
			double rx = x[j] - x[iGlobal];
			double ry = y[j] - y[iGlobal];
			double rz = z[j] - z[iGlobal];
			
			int Ntri = e[Nt_max * j];

			double G_sum = 0.0, F_sum = 0.0;

			if (iGlobal == j) { // SINGULAR elements
				G_sum = G_diag[iGlobal];
			}
			else {
				for (int t = 1; t <= Ntri; t++) { // REGULAR cycle through all j-adjacent triangles

					double Gauss_sum_G = 0.0;
					double Gauss_sum_F = 0.0;

					for (int k = 0; k < NGauss; k++) { // cycle through all Gauss pts of a triangle

						double x_ijk = r_ijk_x[Nt_max * NGauss * j + NGauss * t + k] - x[iGlobal];
						double y_ijk = r_ijk_y[Nt_max * NGauss * j + NGauss * t + k] - y[iGlobal];
						double z_ijk = r_ijk_z[Nt_max * NGauss * j + NGauss * t + k] - z[iGlobal];

						double norm = sqrt(x_ijk * x_ijk + y_ijk * y_ijk + z_ijk * z_ijk);

						Gauss_sum_F += etha_1[k] / (norm * norm * norm) * w[k];
						Gauss_sum_G += etha_1[k] / norm * w[k];
					}

					// triangle area
					double At = A_t[Nt_max * j + t];

					double nx = nt_x[Nt_max * j + t];
					double ny = nt_y[Nt_max * j + t];
					double nz = nt_z[Nt_max * j + t];

					double K_ijt = fabs(rx * nx + ry * ny + rz * nz);

					G_sum += At * Gauss_sum_G;
					F_sum += At * K_ijt * Gauss_sum_F;
				}
				// ------------- end regular triangle cycle ------------------
			}

			double G_ij = 1.0 / (4 * M_PI) * G_sum;
			F_local[i][j] = 1.0 / (4 * M_PI) * F_sum;

			rhs_sum += G_ij * q[j];
		}

		rhs_local[i] = rhs_sum;
	}

	iGlobal = 0;
	for (int i = 0; i < nlocal; i++) { // filling in the diagonal terms of F
		iGlobal = i + istart;
		double F_sum = 0.0;

		for (int j = 0; j < N; j++) {
			if (i == j) continue;

			F_sum += F_local[i][j];
		}

		F_diag[iGlobal] = F_local[i][i] = 1.0 - F_sum;
	}

	MPI_Allgather(rhs_local, nlocal, MPI_DOUBLE, rhs, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);

	auto endMatrixGen = std::chrono::high_resolution_clock::now();

	printArray1FromTo("F_diag", F_diag, istart, istart + nlocal, 4);
	printArray1FromTo("rhs", rhs, istart, istart + nlocal, 4);

	// ============================================================================================================
	// ========== END BOUNDARY ELEM METHOD ========================================================================

	double* u = new double[N]; // potential solution

	// Bi-CGSTAB solve:
	// =========================================================================================
	// ========== Bi-CGSTAB Implementation inside main() =======================================
	auto startBi_CGSTAB = std::chrono::high_resolution_clock::now();

	if (myrank == 0) {
		std::cout << "==================================================" << std::endl;
		std::cout << "----------- Initializing Bi-CGSTAB Method --------" << std::endl;
	}

	// x0 = (1000,1000,...,1000)
	for (int i = 0; i < N; i++) x_curr[i] = 1000.;
	// r0 = b - A x0
	// choose rp0 such that <r0, rp0> != 0
	// p0 = r0

	// local vector packets:
	double* r_curr_local = new double[nlocal];
	double* rp0_local = new double[nlocal];
	double* p_curr_local = new double[nlocal];

	for (int i = 0; i < nlocal; i++) {
		iGlobal = i + istart;
		r_curr_local[i] = rhs[iGlobal];
		for (int j = 0; j < N; j++) {
			r_curr_local[i] -= F_local[i][j] * x_curr[j];
		}
		rp0_local[i] = r_curr_local[i] + 100;
		p_curr_local[i] = r_curr_local[i];
	}

	// gathering vector packets:
	MPI_Allgather(r_curr_local, nlocal, MPI_DOUBLE, r_curr, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(rp0_local, nlocal, MPI_DOUBLE, rp0, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(p_curr_local, nlocal, MPI_DOUBLE, p_curr, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);

	if (myrank == 0) std::cout << "------------ Launching iterations ----------------" << std::endl;
	// begin iterations
	for (int k = 0; k < maxIter; k++) {
		if (myrank == 0) std::cout << "::: iter : " << k << std::endl;

		// alpha[k] = <r[k], rp0> / <Ap[k], rp0>

		double num = 0.; double den = 0.;
		for (int i = 0; i < N; i++) {
			num += r_curr[i] * rp0[i];
		}

		double* tmp_local = new double[nlocal];
		for (int i = 0; i < nlocal; i++) {
			tmp_local[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp_local[i] += F_local[i][j] * p_curr[j];
			}
		}
		MPI_Allgather(tmp_local, nlocal, MPI_DOUBLE, tmp, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);

		for (int i = 0; i < N; i++) {
			den += tmp[i] * rp0[i];
		}

		alpha = num / den;

		// s[k] = r[k] - alpha[k] * A p[k]

		for (int i = 0; i < N; i++) {
			s[i] = r_curr[i] - alpha * tmp[i];
		}

		double norm = vectorNorm(s);

		if (myrank == 0) std::cout << "||s|| = " << norm << std::endl;
		if (norm < tol) {
			// x[k + 1] = x[k] + alpha[k] * p[k]

			for (int i = 0; i < N; i++) {
				x_next[i] = x_curr[i] + alpha * p_curr[i];
			}

			if (myrank == 0) std::cout << "||s|| < tol = " << tol << ", exiting iterations" << std::endl;
			break;
		}

		// omega[k] = <A s[k], s[k]> / <A s[k], A s[k]>

		num = 0; den = 0;
		for (int i = 0; i < nlocal; i++) {
			tmp_local[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp_local[i] += F_local[i][j] * s[j];
			}
		}
		MPI_Allgather(tmp_local, nlocal, MPI_DOUBLE, tmp, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);

		for (int i = 0; i < N; i++) {
			num += tmp[i] * s[i];
			den += tmp[i] * tmp[i];
		}
		omega = num / den;

		// x[k + 1] = x[k] + alpha[k] * p[k] + omega[k] * s[k]
		// r[k + 1] = s[k] - omega[k] * A s[k]

		for (int i = 0; i < N; i++) {
			x_next[i] = x_curr[i] + alpha * p_curr[i] + omega * s[i];
			r_next[i] = s[i] - omega * tmp[i];
		}

		norm = vectorNorm(r_next);
		if (myrank == 0) std::cout << "||r[k + 1]|| = " << norm << std::endl;
		if (norm < tol) {
			if (myrank == 0) std::cout << "||r[k + 1]|| < tol = " << tol << ", exiting iterations" << std::endl;
			delete[] tmp_local;
			break;
		}

		// beta[k] = (alpha[k] / omega[k]) * <r[k + 1], rp0> / <r[k], rp0>

		num = 0; den = 0;
		for (int i = 0; i < N; i++) {
			num += r_next[i] * rp0[i];
			den += r_curr[i] * rp0[i];
		}

		beta = (alpha / omega) * num / den;

		// p[k + 1] = r[k + 1] + beta[k] * (p[k] - omega[k] * A p[k])

		for (int i = 0; i < nlocal; i++) {
			tmp_local[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp_local[i] += F_local[i][j] * p_curr[j];
			}
		}
		MPI_Allgather(tmp_local, nlocal, MPI_DOUBLE, tmp, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);
		delete[] tmp_local;

		for (int i = 0; i < N; i++) {
			p_next[i] = r_next[i] + beta * (p_curr[i] - omega * tmp[i]);
		}

		norm = fabs(vectorDot(r_next, rp0));
		if (myrank == 0) std::cout << "|< r[k + 1], rp0 >| = " << norm << std::endl;
		if (norm < tol) {
			// rp0 = r[k + 1]; p[k + 1] = r[k + 1]

			for (int i = 0; i < N; i++) {
				rp0[i] = r_next[i]; p_next[i] = r_next[i];
			}
		}
		// current = next

		for (int i = 0; i < N; i++) {
			x_curr[i] = x_next[i];
			r_curr[i] = r_next[i];
			p_curr[i] = p_next[i];
		}

		if (myrank == 0) std::cout << "===> finishing iter " << k << std::endl;
	}

	for (int i = 0; i < N; i++) u[i] = x_next[i];

	auto endBi_CGSTAB = std::chrono::high_resolution_clock::now();
	// ========================= Solution done ===============================================

	if (myrank == 0) {
		// print solution
		printArray1("u", u, 4);

		// print exact solution
		printExactSol("uExact", 4);

		std::fstream dataOut;
		dataOut.open("data.dat", std::fstream::out);
		if (!dataOut.is_open()) {
			std::cout << "unable to open output file data.dat!" << std::endl;
			return 0;
		}

		for (int i = 0; i < N; i++) {
			dataOut << B[i] << " " << L[i] << " " << u[i] << std::endl;
		}

		dataOut.close();

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
		std::cout << "..........................................................................." << std::endl;
		std::chrono::duration<double> elapsedMatrixGen = (endMatrixGen - startMatrixGen);
		std::cout << "generating BEM system matrix :    " << elapsedMatrixGen.count() << " s" << std::endl;
		std::cout << ".... Bi-CGSTAB: .........................................................." << std::endl;
		std::chrono::duration<double> elapsedBi_CGSTAB = (endBi_CGSTAB - startBi_CGSTAB);
		std::cout << "Bi-CGSTAB solution :    " << elapsedBi_CGSTAB.count() << " s" << std::endl;
		std::cout << "--------------------------------------------------------------------------" << std::endl;
	}	

	MPI_Finalize();

	// clean up

	delete[] nt_x; delete[] nt_y; delete[] nt_z;
	delete[] G_diag; delete[] A_t;
	delete[] r_ijk_x; delete[] r_ijk_y; delete[] r_ijk_z;

	// delete[] F_diag;
	delete[] r_curr_local; delete[] rp0_local; delete[] p_curr_local;

	delete[] x_curr; delete[] x_next;
	delete[] r_curr; delete[] r_next;
	delete[] p_curr; delete[] p_next;
	delete[] s; delete[] tmp;

	
	delete[] x; delete[] y; delete[] z;
	delete[] q; delete[] d2U;

	for (int i = 0; i < N; i++) delete[] F_local[i];
	delete[] F_local;

	return 1;
}