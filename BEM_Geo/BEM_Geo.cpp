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
double R = 6378000.; // Earth radius [m]
double r_remaining = 500000; // [m]
double GM = 389600.5; // [km^3 * s^(-1)]
double uExact = GM / R; // exact surface potential
double qExact = GM / (R * R); // surface acceleration

const char* dataFilename = "BL-902.dat";
// const char* dataFilename = "BL-1298.dat";
// const char* dataFilename = "BL-3602.dat";
// const char* dataFilename = "BL-8102.dat";
const char* elemDataFilename = "elem_902.dat";
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

void loadPointData(
	double* B, double* L, double* H,
	double* x, double* y, double* z,
	double* sx, double* sy, double* sz,
	double* q, double* d2U
) {
	printf("Loading data ... \n");

	std::fstream dataFile;
	dataFile.open(dataFilename, std::fstream::in);

	if (!dataFile.is_open()) {
		printf("Unable to open file %s\n", dataFilename);
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

			sx[i] = (R - r_remaining) * cos(B[i] * M_PI / 180) * cos(L[i] * M_PI / 180);
			sy[i] = (R - r_remaining) * cos(B[i] * M_PI / 180) * sin(L[i] * M_PI / 180);
			sz[i] = (R - r_remaining) * sin(B[i] * M_PI / 180);

			x[i] = (R + H[i]) * cos(B[i] * M_PI / 180) * cos(L[i] * M_PI / 180);
			y[i] = (R + H[i]) * cos(B[i] * M_PI / 180) * sin(L[i] * M_PI / 180);
			z[i] = (R + H[i]) * sin(B[i] * M_PI / 180);

			q[i] = 0.00001 * Q;
			d2U[i] = D2U;

			// std::cout << B << " " << L << " " << H << " " << Q << " " << D2U << std::endl;
			// std::cout << sx[i] << " " << sy[i] << " " << sz[i] << " " << x[i] << " " << y[i] << " " << z[i] << " " << q[i] << " " << d2U[i] << std::endl << std::endl;
			i++;
		}

		dataFile.close();
	}
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

void getMatrices(double** dG, double** G, double* x, double* y, double* z, double* sx, double* sy, double* sz) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double dx = x[i] - sx[j];
			double dy = y[i] - sy[j];
			double dz = z[i] - sz[j];

			double d_norm = sqrt(dx * dx + dy * dy + dz * dz);

			double nx = sx[i] / (R - r_remaining);
			double ny = sy[i] / (R - r_remaining);
			double nz = sz[i] / (R - r_remaining);

			double dot = dx * nx + dy * ny + dz * nz;

			dG[i][j] = dot / (4 * M_PI * d_norm * d_norm * d_norm); // dG[i][j]/dn[i]
			G[i][j] = 1 / (4 * M_PI * d_norm);
			if (isnan(dG[i][j])) {
				std::cout << "NAN!: dG[" << i << "][" << j << "] : d_norm = " << d_norm << std::endl;
				std::cout << "dx = " << dx << ", dy = " << dy << ", dz = " << dz << std::endl;
			}
			if (isnan(G[i][j])) {
				std::cout << "NAN!: G[" << i << "][" << j << "] : d_norm = " << d_norm << std::endl;
				std::cout << "dx = " << dx << ", dy = " << dy << ", dz = " << dz << std::endl;
			}
		}
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

void mmultVector(double** A, double* x, double* result) {
	for (int i = 0; i < N; i++) {
		result[i] = 0.;
		for (int j = 0; j < N; j++) {
			result[i] += A[i][j] * x[j];
		}
	}
}

void smultVector(double s, double* x, double* result) {
	for (int i = 0; i < N; i++) result[i] = x[i] * s;
}

void subVectors(double* a, double* b, double* result) {
	for (int i = 0; i < N; i++) result[i] = a[i] - b[i];
}

void addVectors(double* a, double* b, double* result) {
	for (int i = 0; i < N; i++) result[i] = a[i] + b[i];
}

void copyVector(double* original, double* target) {
	for (int i = 0; i < N; i++) target[i] = original[i];
}

void Bi_CGSTAB_solve(double** A, double* b, double* x) {
	// ctrl. constants
	int maxIter = 1000;
	double tol = 1e-5;

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
	double alpha, beta, omega;

	// x0 = (1,1,...,1)
	for (int i = 0; i < N; i++) x_curr[i] = 1000.;
	// r0 = b - A x0
	// choose rp0 such that <r0, rp0> != 0
	// p0 = r0
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
		std::cout << "::: iter : " << k << std::endl;
		// alpha[k] = <r[k], rp0> / <Ap[k], rp0>
		mmultVector(A, p_curr, tmp);
		alpha = vectorDot(r_curr, rp0) / vectorDot(tmp, rp0);
		// s[k] = r[k] - alpha[k] * A p[k]
		smultVector(alpha, tmp, tmp);
		subVectors(r_curr, tmp, s);

		std::cout << "||s|| = " << vectorNorm(s) << std::endl;
		if (vectorNorm(s) < tol) {
			// x[k + 1] = x[k] + alpha[k] * p[k]
			smultVector(alpha, p_curr, tmp);
			addVectors(x_curr, tmp, x_next);
			std::cout << "||s|| < tol = " << tol << ", exiting iterations" << std::endl;
			break;
		}

		// omega[k] = <A s[k], s[k]> / <A s[k], A s[k]>
		mmultVector(A, s, tmp);
		omega = vectorDot(tmp, s) / vectorDot(tmp, tmp);
		// x[k + 1] = x[k] + alpha[k] * p[k] + omega[k] * s[k]
		smultVector(alpha, p_curr, tmp);
		addVectors(x_curr, tmp, tmp);
		smultVector(omega, s, tmp1);
		addVectors(tmp, tmp1, x_next);
		// printArray1("x", x_next, 5);
		// r[k + 1] = s[k] - omega[k] * A s[k]
		mmultVector(A, s, tmp);
		smultVector(omega, tmp, tmp);
		subVectors(s, tmp, r_next);
		std::cout << "||r[k + 1]|| = " << vectorNorm(r_next) << std::endl;
		if (vectorNorm(r_next) < tol) {
			std::cout << "||r[k + 1]|| < tol = " << tol << ", exiting iterations" << std::endl;
			break;
		}
		// beta[k] = (alpha[k] / omega[k]) * <r[k + 1], rp0> / <r[k], rp0>
		beta = (alpha / omega) * vectorDot(r_next, rp0) / vectorDot(r_curr, rp0);
		// p[k + 1] = r[k + 1] + beta[k] * (p[k] - omega[k] * A p[k])
		mmultVector(A, p_curr, tmp);
		smultVector(omega, tmp, tmp);
		subVectors(p_curr, tmp, tmp);
		smultVector(beta, tmp, tmp);
		addVectors(r_next, tmp, p_next);
		std::cout << "|< r[k + 1], rp0 >| = " << fabs(vectorDot(r_next, rp0)) << std::endl;
		if (fabs(vectorDot(r_next, rp0)) < tol) {
			// rp0 = r[k + 1]; p[k + 1] = r[k + 1]
			// std::cout << "|< r[k + 1], rp0 >| < tol = " << tol << ", copying r[k + 1] to rp0 and p[k + 1]" << std::endl;
			copyVector(r_next, rp0);
			copyVector(r_next, p_next);
		}
		// current = next
		copyVector(x_next, x_curr);
		copyVector(r_next, r_curr);
		copyVector(p_next, p_curr);
		std::cout << "===> finishing iter " << k << std::endl;
	}

	copyVector(x_next, x); // result: x = x_next

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

int main() {
	// point data

	double* B = new double[N];
	double* L = new double[N];
	double* H = new double[N];

	double* x = new double[N];
	double* y = new double[N];
	double* z = new double[N];

	double* sx = new double[N];
	double* sy = new double[N];
	double* sz = new double[N];

	double* q = new double[N];
	double* d2U = new double[N];

	// elem data

	int* e = new int[(size_t)Nt_max * N];

	auto startLoad = std::chrono::high_resolution_clock::now();
	loadPointData(B, L, H, x, y, z, sx, sy, sz, q, d2U);
	loadElemData(e);
	auto endLoad = std::chrono::high_resolution_clock::now();

	auto startPrint3 = std::chrono::high_resolution_clock::now();
	printArrayVector3("x", x, y, z, 5);
	auto endPrint3 = std::chrono::high_resolution_clock::now();
	printArrayVector3("s", sx, sy, sz, 5);

	auto startPrint1 = std::chrono::high_resolution_clock::now();
	printArray1("q", q, 5);
	auto endPrint1 = std::chrono::high_resolution_clock::now();
	printArray1("d2U/dn2", d2U, 5);

	printVertNeighborArray("e", e, 7);



	// ==============================================================================================================
	// ========== BOUNDARY ELEMENT METHOD ===========================================================================


	// Generating system matrices from point and elem data:

	auto startMatrixGen = std::chrono::high_resolution_clock::now();

	double** F = new double* [N];
	double** G = new double* [N];
	for (int i = 0; i < N; i++) {
		F[i] = new double[N];
		G[i] = new double[N];
	}

	int Ntri, i, j, t, k;
	double rx, ry, rz; // r_ij distance coords
	double nx, ny, nz;
	double x_ijk, y_ijk, z_ijk; // Gauss pts
	double e0x, e0y, e0z, e1x, e1y, e1z, e2x, e2y, e2z; // triangle edges
	double alpha_t, beta_t, l0_t, l1_t, l2_t; // singular triangle features (angle at vertex,  angle at adjacent vertex,  opposing edge length)
	double A_t, K_ijt, norm;
	double G_sum, F_sum, Gauss_sum;
	double t0x, t0y, t0z, t1x, t1y, t1z, t2x, t2y, t2z;
	int j1, j2;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			rx = x[j] - x[i];
			ry = y[j] - y[i];
			rz = z[j] - z[i];
			
			Ntri = e[Nt_max * j];

			G_sum = 0.0; F_sum = 0.0;

			if (i == j) {
				for (t = 1; t <= Ntri; t++) { // SINGULAR cycle through all j=i-adjacent triangles

					// adjacent vertex indices
					j1 = e[Nt_max * j + t];
					j2 = e[Nt_max * j + t % Ntri + 1];

					// compute edge coords
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
					
					alpha_t = acos(-(e0x * e2x + e0y * e2y + e0z * e2z) / (l0_t * l2_t));
					beta_t = acos(-(e0x * e1x + e0y * e1y + e0z * e1z) / (l0_t * l1_t));

					// tri verts
					t0x = x[j];		t0y = y[j];		t0z = z[j];
					t1x = x[j1];	t1y = y[j1];	t1z = z[j1];
					t2x = x[j2];	t2y = y[j2];	t2z = z[j2];

					// tri normal
					nx = (t1y - t0y) * (t2z - t0z) - (t1z - t0z) * (t2y - t0y);
					ny = (t1z - t0z) * (t2x - t0x) - (t1x - t0x) * (t2z - t0z);
					nz = (t1x - t0x) * (t2y - t0y) - (t1y - t0y) * (t2x - t0x);

					norm = sqrt(nx * nx + ny * ny + nz * nz);
					A_t = 0.5 * norm; // triangle area

					assert(norm != 0.0);

					G_sum += A_t / l1_t * log(tan(0.5 * (alpha_t + beta_t)) / tan(0.5 * beta_t));
				}
			}
			else {
				for (t = 1; t <= Ntri; t++) { // REGULAR cycle through all j-adjacent triangles

					// adjacent vertex indices
					j1 = e[Nt_max * j + t];
					j2 = e[Nt_max * j + t % Ntri + 1];
					Gauss_sum = 0.0;

					for (k = 0; k < NGauss; k++) { // cycle through all Gauss pts of a triangle

						x_ijk = etha_1[k] * x[j] + etha_2[k] * x[j1] + etha_3[k] * x[j2] - x[i];
						y_ijk = etha_1[k] * y[j] + etha_2[k] * y[j1] + etha_3[k] * y[j2] - y[i];
						z_ijk = etha_1[k] * z[j] + etha_2[k] * z[j1] + etha_3[k] * z[j2] - z[i];

						norm = sqrt(x_ijk * x_ijk + y_ijk * y_ijk + z_ijk * z_ijk);

						assert(norm != 0.0);

						Gauss_sum += etha_1[k] / (norm * norm * norm) * w[k];
					}

					// tri verts
					t0x = x[j];		t0y = y[j];		t0z = z[j];
					t1x = x[j1];	t1y = y[j1];	t1z = z[j1];
					t2x = x[j2];	t2y = y[j2];	t2z = z[j2];

					// tri normal
					nx = (t1y - t0y) * (t2z - t0z) - (t1z - t0z) * (t2y - t0y);
					ny = (t1z - t0z) * (t2x - t0x) - (t1x - t0x) * (t2z - t0z);
					nz = (t1x - t0x) * (t2y - t0y) - (t1y - t0y) * (t2x - t0x);

					norm = sqrt(nx * nx + ny * ny + nz * nz);
					A_t = 0.5 * norm; // triangle area

					assert(norm != 0.0);

					nx /= norm;
					ny /= norm;
					nz /= norm;

					K_ijt = fabs(rx * nx + ry * ny + rz * nz);

					G_sum += A_t * Gauss_sum;
					F_sum += A_t * K_ijt * Gauss_sum;

					assert(F_sum >= 0.0);
				}
				// ------- end regular triangle cycle --------
			}

			G[i][j] = 1.0 / (4 * M_PI) * G_sum;
			F[i][j] = 1.0 / (4 * M_PI) * F_sum;
		}
	}

	for (i = 0; i < N; i++) { // filling in the diagonal terms of F

		F_sum = 0.0;

		for (j = 0; j < N; j++) {
			if (i == j) continue;

			F_sum += F[i][j];
		}

		F[i][i] = 1.0 - F_sum;
	}

	auto endMatrixGen = std::chrono::high_resolution_clock::now();

	// ================ rhs = G . q ===============================================

	double* rhs = new double[N];

	for (i = 0; i < N; i++) {
		G_sum = 0.0;
		for (j = 0; j < N; j++) {
			G_sum += G[i][j] * q[j];
		}
		rhs[i] = G_sum;
	}


	auto startMatrixPrint = std::chrono::high_resolution_clock::now();
	printArray2("F", F, 4);
	auto endMatrixPrint = std::chrono::high_resolution_clock::now();
	printArray2("G", G, 4);

	// ============================================================================================================
	// ========== END BOUNDARY ELEM METHOD ========================================================================

	double* u = new double[N]; // potential solution

	// Bi-CGSTAB solve:
	// NOTE: This has to be done inside the main() function
	// auto startBi_CGSTAB = std::chrono::high_resolution_clock::now();
	// Bi_CGSTAB_solve(dG, q, alphas); // you had a good life, but I can't use you! :'(
	// auto endBi_CGSTAB = std::chrono::high_resolution_clock::now();

	// =========================================================================================
	// ========== Bi-CGSTAB Implementation inside main() =======================================
	auto startBi_CGSTAB = std::chrono::high_resolution_clock::now();

	// ctrl. constants
	int maxIter = 100;
	double tol = 1e-5;

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
	double alpha, beta, omega;

	// x0 = (1000,1000,...,1000)
	for (int i = 0; i < N; i++) x_curr[i] = 1000.;
	// r0 = b - A x0
	// choose rp0 such that <r0, rp0> != 0
	// p0 = r0
	for (int i = 0; i < N; i++) {
		r_curr[i] = q[i];
		for (int j = 0; j < N; j++) {
			r_curr[i] -= F[i][j] * x_curr[j];
		}
		rp0[i] = r_curr[i] + 100;
		p_curr[i] = r_curr[i];
	}
	std::cout << "==================================================" << std::endl;
	std::cout << "----------- Initializing Bi-CGSTAB Method --------" << std::endl;
	printArray2("systemMatrix", F, 4);
	printArray1("systemRhs", q, 5);
	printArray1("x0", x_curr, 2);
	printArray1("r0", r_curr, 5);

	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "------------ Launching iterations ----------------" << std::endl;
	// begin iterations
	for (int k = 0; k < maxIter; k++) {
		std::cout << "::: iter : " << k << std::endl;

		// alpha[k] = <r[k], rp0> / <Ap[k], rp0>

		double num = 0.; double den = 0.;
		for (int i = 0; i < N; i++) {
			tmp[i] = 0.;
			for (int j = 0; j < N; j++) {
				tmp[i] += F[i][j] * p_curr[j];
			}
			num += r_curr[i] * rp0[i];
			den += tmp[i] * rp0[i];
		}
		alpha = num / den;

		// s[k] = r[k] - alpha[k] * A p[k]

		for (int i = 0; i < N; i++) {
			s[i] = r_curr[i] - alpha * tmp[i];
		}

		double norm = vectorNorm(s);
		std::cout << "||s|| = " << norm << std::endl;
		if (norm < tol) {
			// x[k + 1] = x[k] + alpha[k] * p[k]

			for (int i = 0; i < N; i++) {
				x_next[i] = x_curr[i] + alpha * p_curr[i];
			}

			std::cout << "||s|| < tol = " << tol << ", exiting iterations" << std::endl;
			break;
		}

		// omega[k] = <A s[k], s[k]> / <A s[k], A s[k]>

		num = 0; den = 0;
		for (int i = 0; i < N; i++) {
			tmp[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp[i] += F[i][j] * s[j];
			}
			num += tmp[i] * s[i];
			den += tmp[i] * tmp[i];
		}
		omega = num / den;

		// x[k + 1] = x[k] + alpha[k] * p[k] + omega[k] * s[k]

		for (int i = 0; i < N; i++) {
			x_next[i] = x_curr[i] + alpha * p_curr[i] + omega * s[i];
		}

		// r[k + 1] = s[k] - omega[k] * A s[k]

		for (int i = 0; i < N; i++) {
			r_next[i] = s[i] - omega * tmp[i];
		}

		norm = vectorNorm(r_next);
		std::cout << "||r[k + 1]|| = " << norm << std::endl;
		if (norm < tol) {
			std::cout << "||r[k + 1]|| < tol = " << tol << ", exiting iterations" << std::endl;
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

		for (int i = 0; i < N; i++) {
			tmp[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp[i] += F[i][j] * p_curr[j];
			}
			p_next[i] = r_next[i] + beta * (p_curr[i] - omega * tmp[i]);
		}

		norm = fabs(vectorDot(r_next, rp0));
		std::cout << "|< r[k + 1], rp0 >| = " << norm << std::endl;
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

		std::cout << "===> finishing iter " << k << std::endl;
	}

	// result: x = x_next
	for (int i = 0; i < N; i++) u[i] = x_next[i];

	// clean up
	delete[] x_curr; delete[] x_next;
	delete[] r_curr; delete[] r_next;
	delete[] p_curr; delete[] p_next;
	delete[] s; delete[] tmp; delete[] tmp1;

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
	std::chrono::duration<double> elapsedPrint3 = (endPrint3 - startPrint3);
	std::cout << "printing vector 3 array :    " << elapsedPrint3.count() << " s" << std::endl;
	std::chrono::duration<double> elapsedPrint1 = (endPrint1 - startPrint1);
	std::cout << "printing vector 1 array :    " << elapsedPrint3.count() << " s" << std::endl;
	std::chrono::duration<double> elapsedPrintMatrix = (endMatrixPrint - startMatrixPrint);
	std::cout << "printing matrix array :    " << elapsedPrintMatrix.count() << " s" << std::endl;
	std::cout << "..........................................................................." << std::endl;
	std::chrono::duration<double> elapsedMatrixGen = (endMatrixGen - startMatrixGen);
	std::cout << "generating system matrix :    " << elapsedMatrixGen.count() << " s" << std::endl;
	std::cout << ".... Bi-CGSTAB: .........................................................." << std::endl;
	std::chrono::duration<double> elapsedBi_CGSTAB = (endBi_CGSTAB - startBi_CGSTAB);
	std::cout << "Bi-CGSTAB solution :    " << elapsedBi_CGSTAB.count() << " s" << std::endl;
	std::cout << "--------------------------------------------------------------------------" << std::endl;


	// clean up
	delete[] x; delete[] y; delete[] z;
	delete[] sx; delete[] sy; delete[] sz;
	delete[] q; delete[] d2U;

	for (int i = 0; i < N; i++) {
		delete[] F[i]; delete[] G[i];
	}
	delete[] F; delete[] G;

	return 1;
}