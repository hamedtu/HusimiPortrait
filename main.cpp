/*
This program computes Floquetstates for the periodically driven Bose-Hubbard Dimer
using a "time-slicing" method. The Hamiltonian is approximated as constant within
each small time-slice (∆t), and the evolution operator for each slice is
calculated as exp(-iH∆t). The total evolution operator is obtained by multiplying
the operators for each time slice. The method requires diagonalizing the Hamiltonian
at each time step to obtain eigenvalues and eigenstates. The computation is parallelized
using MPI for efficient handling of large systems with many time slices. This method
is most effective when the driving Hamiltonian is much smaller than the undriven part.

Program Arguments:
param[0] - System parameter (Nk/Omega)
param[1] - Scaled driving strength
param[2] - Scaled driving frequency (small omega)
N        - Number of particles
L        - Number of time slices

Author: Steffen Seligmann
Paper: "Degree of Simplicity of Floquet States of a Periodically Driven Bose-Hubbard Dimer"
Supplementary code.
*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <Eigen/Dense>
#include <iomanip>
#include <sstream>
#include <string>
#include <mpi.h>

using namespace std;
using namespace Eigen;

// Complex number definition
complex<double> i(0, 1); // Imaginary unit for complex number operations
const double PI = 3.141592653589793238463; // Define the value of Pi

// Type definitions for convenience
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RVectorXd; // Real vector of dynamic size
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> RMatrixXd; // Real matrix of dynamic size
typedef Eigen::Matrix<complex<double>, 1, Eigen::Dynamic> CVectorXd; // Complex vector of dynamic size
typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> CMatrixXd; // Complex matrix of dynamic size

// Forward declaration of the spectral function, which is used later in the code.
CMatrixXd spectral(RVectorXd E_n, CMatrixXd V_n, double deltaT, int N);

int main(int argc, char *argv[])
{
    // Initialize MPI environment
    MPI_Init(NULL, NULL); // Initializes the MPI environment
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Get the total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // Get the rank of the current process

    // Parse parameters from command-line arguments
    double param[3] = {stod(argv[1]), stod(argv[2]), stod(argv[3])}; // Store the first three parameters
    const int N = stoi(argv[4]); // Size of the system (e.g., number of states)
    int L = stoi(argv[5]); // Number of time steps or iterations
    double tLim[2] = {0.0, (PI / (2 * param[2]))}; // Define time limits based on parameter [2]

    // Generate time values for simulation using the LinSpaced function from Eigen
    RVectorXd t = RVectorXd::LinSpaced(L, 1, 0); // Linearly spaced time points

    // Transform time points based on parameter [2]
    for (int l = 0; l < L; l++)
    {
        t(l) = acos(t(l)) / param[2]; // Adjusting time based on a physical model
    }

    // Hamiltonian matrix without the drive (initial system)
    RMatrixXd H_0 = RMatrixXd::Zero(N + 1, N + 1); // Initialize zero matrix of size (N+1)x(N+1)
    H_0.diagonal() = param[0] * (RVectorXd::LinSpaced(N + 1, -N, N).array().pow(2) / (2 * N) + (N - 2) / 2.0); // Diagonal elements based on energy

    // Fill off-diagonal elements of the Hamiltonian matrix
    H_0.block(1, 0, N, N).diagonal() = -1 * (RVectorXd::LinSpaced(N, 1, N).array() * RVectorXd::LinSpaced(N, N, 1).array()).sqrt() / 2;
    H_0.block(0, 1, N, N).diagonal() = -1 * (RVectorXd::LinSpaced(N, 1, N).array() * RVectorXd::LinSpaced(N, N, 1).array()).sqrt() / 2;

    // Drive coefficient matrix (interacts with H_0 to simulate time evolution)
    RMatrixXd H_1 = (RVectorXd::LinSpaced(N + 1, -N, N) * param[1]).asDiagonal(); // Diagonal matrix based on the second parameter (param[1])
    CMatrixXd U_T = CMatrixXd::Identity(N + 1, N + 1); // Identity matrix (initial state)
    Eigen::SelfAdjointEigenSolver<RMatrixXd> H_SOL; // Solver for the Hamiltonian eigenvalues and eigenvectors

    // MPI parallelization setup: divide work among processes
    int chunk_size = floor((L-1) / world_size); // Calculate the chunk size for each process
    int extra_iterations = (L -1) % world_size; // Handle any leftover iterations
    const int localIterations = chunk_size + (world_rank < (L-1) % world_size ? 1 : 0); // Calculate local iterations for each process
    int start_index = world_rank * chunk_size + std::min(world_rank, extra_iterations); // Starting index for the current process
    int end_index = start_index + localIterations; // Ending index for the current process

    // Start measuring time for the computation
    auto start = std::chrono::steady_clock::now();

    // Main time-stepping loop: parallelized across MPI processes
    for (int i = start_index; i < end_index; ++i)
    {
        // Compute the momentary Hamiltonian (with time dependence)
        H_SOL.compute(H_0 + H_1 * cos(param[2] * t(i))); // Diagonalize the Hamiltonian
        U_T = spectral(H_SOL.eigenvalues(), H_SOL.eigenvectors(), t(i + 1) - t(i), N) * U_T; // Update the evolution matrix
    }

    const int matrix_size = (N + 1) * (N + 1); // Total size of the matrix (for MPI communication)
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processes before sending/receiving data

    if (world_rank == 0)
    {
        // Root process: receive matrices from other processes and multiply them
        for (int i = 1; i < world_size; ++i)
        {
            // Receive matrix data from process i
            CMatrixXd received_matrix(N + 1, N + 1);
            MPI_Recv(received_matrix.data(), matrix_size * sizeof(std::complex<double>), MPI_BYTE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // Multiply received matrix with the accumulated result
            U_T = received_matrix * U_T;
        }
    }
    else
    {
        // Non-root processes: send the matrix data to the root process
        MPI_Send(U_T.data(), matrix_size * sizeof(std::complex<double>), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD); // Ensure synchronization again before proceeding

    if (world_rank == 0)
    {
        // Root process: Use symmetry to simplify the result (e.g., transpose and multiply matrices)
        U_T = U_T(seqN(last, N + 1, fix<-1>), seqN(last, N + 1, fix<-1>)).transpose() * U_T;
        U_T = U_T.transpose() * U_T; // Symmetry operation on the evolution matrix

        // Eigenvalue decomposition of the final matrix
        ComplexEigenSolver<CMatrixXd> U_SOL;
        U_SOL.compute(U_T); // Compute the eigenvectors and eigenvalues of the final matrix

        // Measure elapsed time for the computation
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

        // Convert numerical parameters to strings with specific precision for file naming
        std::stringstream ssAlpha, ssMu, ssW, ssTolerance;
        ssAlpha << std::fixed << std::setprecision(2) << param[0];
        ssMu << std::fixed << std::setprecision(2) << param[1];
        ssW << std::fixed << std::setprecision(2) << param[2];

        std::string strAlpha = ssAlpha.str();
        std::string strMu = ssMu.str();
        std::string strW = ssW.str();

        // Construct file names based on parameters
        string file_real = "TS_eigenvectors_" + strAlpha + "_" + strMu + "_" + strW + "_" + to_string(N) + "_" + to_string(L) + "_real.dat";
        string file_imag = "TS_eigenvectors_" + strAlpha + "_" + strMu + "_" + strW + "_" + to_string(N) + "_" + to_string(L) + "_imag.dat";

        // Save real and imaginary parts of eigenvectors to separate files
        ofstream file;
        file.open(file_real);
        file << U_SOL.eigenvectors().real() << endl;
        file.close();

        file.open(file_imag);
        file << U_SOL.eigenvectors().imag() << endl;
        file.close();
    }

    // Finalize the MPI environment
    MPI_Finalize();

    return 0; // Return success
}

/**
  Spectral decomposition function

  This function calculates the time evolution operator matrix for the system at a specific time step,
  based on the spectral decomposition, meaning eigenvalues and eigenvectors of the Hamiltonian at time t.

  E_n    -  Eigenvalues of the Hamiltonian.
  V_n    -  Eigenvectors of the Hamiltonian.
  deltaT -  Time step for the evolution.
  N      -  Size of the system.
 */
CMatrixXd spectral(RVectorXd E_n, CMatrixXd V_n, double deltaT, int N)
{
    // Initialize the spectral matrix
    CMatrixXd spectralM = CMatrixXd::Zero(N + 1, N + 1);

    // Construct the spectral matrix by summing over the eigenvectors and eigenvalues
    for (int k = 0; k < N + 1; k++)
    {
        spectralM.noalias() += (V_n.col(k) * V_n.col(k).conjugate().transpose()) * exp(-i * E_n(k) * deltaT); // Exponentially decay each state based on its energy
    }

    return spectralM; // Return the computed spectral matrix
}