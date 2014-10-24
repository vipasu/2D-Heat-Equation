/*
 * Physics 113 Final Project
 * Fall 2013
 * Vincent Su
 *
 * Implementation of various time-stepping schemes for solving the
 * heat equation in one and two dimensions.
 *
 */

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <string.h>
using namespace std;

#define PI M_PI

void Simulate1D(int nx, int nt, double dt, double total_t,
        double (*phi)(double), double alpha2, double x0, double xL,
        void (*advance)(double **, int, double, int, double, double),
        char* mode);
void Simulate2D(int ns, int nt, double dt, double total_t,
                double (*phi)(double, double), double alpha2,
                double x0, double xL, double y0, double yL,
                void (*advance)(double ***, int, double, int,
                    double, double, double, double), char *mode);
double *exact_1D(double t, int nx, double alpha2);
double **exact_2D(double t, int ns, double alpha2);

void FE_1D_Advance(double **up, int nx, double r, int nsteps, double x0, double xL);
void FE_Matrix_Advance(double **up, int nx, double r, int nsteps, double x0, double xL);
void BE_Matrix_Advance(double **up, int nx, double r, int nsteps, double x0, double xL);
void CN_Matrix_Advance(double **up, int nx, double r, int nsteps, double x0, double xL);

void FE_2D_Advance_Simple(double ***mp, int ns, double r, int nsteps,
        double x0, double xL, double y0, double yL);
void FE_2D_Advance(double ***mp, int ns, double r, int nsteps,
        double x0, double xL, double y0, double yL);
void BE_2D_Advance(double ***mp, int ns, double r, int nsteps,
        double x0, double xL, double y0, double yL);
void CN_2D_Advance(double ***mp, int ns, double r, int nsteps,
        double x0, double xL, double y0, double yL);
void CN_2D_ADI_Advance(double ***mp, int ns, double r, int nsteps,
        double x0, double xL, double y0, double yL);

void compare1DMethods(int nx, int nslices, double dt, double t);

/* Helper functions */
double **matrixMultiply(double **S, double **T, int size);
void matrixAdd(double ***Ap, double **B, int n, int m);
void matrixVectorMultiply(double **T, double **Vp, int n, int m);
void solveTriDiag(double **T, double **up, int n);
double **createZeroMatrix(int n, int m);
double **createImplicitLaplacian(int dim, double r);
double **createExplicitLaplacian(int dim, double r);
double **transposeMatrix(double **A, int n, int m);
void freeMatrix(double **arr, int n, int m);
void printMatrix(double **A, int n, int m);
void printMatrix(double **A, int n, int m, int nd);

double *vecApply(double (*f)(double), double *arr, int len);
double *linspace(double start, double stop, int n);
void vectorAdd(double **up, double *v, int len);
void vectorScale(double **up, double c, int len);
void printVec(double *v, int n, int nd);
void printVec(double *v, int n);

double L2Norm(double *v, double *exact, int dim);
double simple_square(double x);
void validateRoutines();


/* Analytical solutions in 1D and 2D for calibration purposes */
double sinx(double x) {
    return sin (x * PI);
}
double *exact_1D(double t, int nx, double alpha2)
{
    double *x = linspace(0,1,nx);
    double damp_factor = exp(-alpha2 * PI * PI * t);
    double *u_0 = vecApply(sinx, x, nx);
    vectorScale(&u_0, damp_factor, nx);
    free(x);
    return u_0;
}

double sinxsiny(double x, double y) {
    return sin( x * PI) * sin (y * PI);
}

double **exact_2D(double t, int ns, double alpha2)
{
    double **grid = createZeroMatrix(ns, ns);
    double damp_factor = exp(-alpha2 * PI * PI * t/2);
    double ds = 1.0/(ns-1);
    for (int i = 1; i < ns-1; i++) {
        for (int j = 1; j < ns-1; j++) {
            grid[i][j] = damp_factor * sinxsiny(i *ds, j * ds);
        }
    }
    return grid;
}

/* Sample 1D IC's */
double flat(double x) {
    return 1;
}
double linear(double x) {
    return 10 * x;
}
double delta(double x) {
    if (abs(x-.5) <= .0001)
        return 100;
    else
        return 0;
}

/* Sample 2D IC's */
double hotrod(double x, double y) {
    if (fabs(y-.5) < .1)
        return 100;
    else
        return 0;
}
double hotcircle(double x, double y) {
    double width = x-.5,
           height = y-.5;
    if (sqrt(width * width + height * height) < .5)
        return 100;
    else
        return 0;
}
double paraboloid(double x, double y){
    double width = x-.5,
           height = y-.5;
    return 1000 * (width*width + height * height);
}

int main()
{
    int ns = 26, nslice = 10;
    double dt = .00001, t = .1;
    double alpha2 = 1, x0 = 0, xL = 0, y0 = 0, yL = 0;
    char display[] = "display",
         calcError[] = "error";

    /* Example usage */

    Simulate2D(ns, nslice, dt, t, paraboloid, alpha2, x0, xL, y0, yL,
            CN_2D_ADI_Advance, display);


    /* Comparison of errors for different time steps */
    /*
    for (int i = 1; i <= 11; i++ ) {
        compare1DMethods(ns, nslice, dt * pow(2,i), t);
    }
    */

    return 0;
}

/* Function: compare1DMethods
 *
 * Inputs:  nx - # of points on the rod
 *          nslices - # of snapshots to take
 *          dt - time step size
 *          t  - total length of simulation
 *
 *  compare1DMethods serves as a wrapper for running all four
 *  implemententations of simulating the rod and prints out their
 *  errors accordingly. For the purposes of this project, error has
 *  been defined as the sum of the L2 norms at each of the
 *  snapshots.
 */
void compare1DMethods(int nx, int nslices, double dt, double t)
{
    double alpha2 = 1, x0 = 0, xL = 0;
    char mode[] = "error";
    Simulate1D(nx, nslices, dt, t, sinx, alpha2, x0, xL, FE_1D_Advance, mode);
    Simulate1D(nx, nslices, dt, t, sinx, alpha2, x0, xL, FE_Matrix_Advance, mode);
    Simulate1D(nx, nslices, dt, t, sinx, alpha2, x0, xL, BE_Matrix_Advance, mode);
    Simulate1D(nx, nslices, dt, t, sinx, alpha2, x0, xL, CN_Matrix_Advance, mode);
}

/* Function: Simulate1D
 *
 * Inputs:
 *      nx  - # of points to put on the rod
 *      nt  - # of times to either print the temperature or evaluate error,
 *              (slices or snapshots)
 *      dt  - length of time step
 *      total_t - length of time to run the simulation
 *      phi - initial temperature profile
 *      alpha2  - diffusivity constant
 *      x0, xL  - Dirichlet boundary conditions
 *      advance - time-stepping scheme (Forward Euler, etc.)
 *      mode    - whether to display the temperature or evaluate error.
 *
 *  Simulate1D serves as a generic wrapper for modeling the time evolution of
 *  the rod. It initializes the simulation with the appropriate parameters,
 *  handles any mismatches with timing, and calls the desired time-stepping
 *  method.
 */
void Simulate1D(int nx, int nt, double dt, double total_t,
        double (*phi)(double), double alpha2, double x0, double xL,
        void (*advance)(double **, int, double, int, double, double),
        char* mode)
{
    double error = 0;
    double *u_exact,
           *xvals = linspace(0, 1, nx);
    double dx = xvals[1] - xvals[0];
    double r = alpha2 * dt/dx/dx;

    /* Enforce IC's and BC's */
    double *u = vecApply(phi, xvals, nx);
    u[0] = x0; u[nx-1] = xL;

    /* Calculate how much time goes between snapshots and rescale the time steps
     * to fit.
     */
    double timeslice = total_t/nt;
    int nsteps = (int) ceil((timeslice/dt));
    dt = timeslice / nsteps;

    for (int i = 0; i <= nt; i++) {
        if (strcmp(mode, "display") == 0){
            printVec(u, nx);
        } else if (strcmp(mode, "error") == 0){
            u_exact = exact_1D(timeslice * i, nx, alpha2);
            error += L2Norm(u, u_exact, nx);
            free(u_exact);
        }
        advance(&u, nx, r, nsteps, x0, xL);
    }

    if (strcmp(mode, "error") == 0){
        printf("%.16e\n", error);
    }
    printf("\n");
    free(xvals);
    free(u);
}

/* Function: Simulate2D
 *
 * Inputs:
 *      ns  - # of points to put on each dimension of the grid
 *      nt  - # of times to either print the temperature or evaluate error,
 *              (slices or snapshots)
 *      dt  - length of time step
 *      total_t - length of time to run the simulation
 *      phi - initial temperature profile
 *      alpha2  - diffusivity constant
 *      x0, xL, y0, yL  - Dirichlet boundary conditions
 *      advance - time-stepping scheme (forward euler, etc.)
 *      mode    - whether to display the temperature or evaluate error.
 *
 *  Simulate2D has the same functionality as Simulate1D, except that it models
 *  the unit square rather than the unit line.
 */
void Simulate2D(int ns, int nt, double dt, double total_t,
                double (*phi)(double, double), double alpha2,
                double x0, double xL, double y0, double yL,
                void (*advance)(double ***, int, double, int,
                    double, double, double, double), char *mode)
{
    double error = 0,
           ds = 1./(ns-1),
           r = alpha2 * dt/ds/ds;
    double **u_exact,
           **mat = createZeroMatrix(ns, ns);
    for (int i = 1; i < ns-1; i++) {
        for (int j = 1; j < ns-1; j++) {
            mat[i][j] = phi(i* ds,j * ds);
        }
    }
    for (int i = 0; i < ns; i++) {
        mat[i][0] = y0;
        mat[i][ns-1] = yL;
        mat[0][i] = x0;
        mat[ns-1][i] = xL;
    }

    double timeslice = total_t/nt;
    int nsteps = (int) ceil((timeslice/dt));
    dt = timeslice / nsteps;

    for (int i = 0; i <= nt; i++) {
        if (strcmp(mode, "display") == 0){
            printMatrix(mat, ns, ns);
        } else if (strcmp(mode, "error") == 0){
            u_exact = exact_2D(timeslice * i, ns, alpha2);
            for (int i = 0; i < ns; i++) {
                error += L2Norm(mat[i], u_exact[i], ns);
            }
            freeMatrix(u_exact, ns, ns);
        }

        advance(&mat, ns, r, nsteps, x0, xL, y0, yL);
    }

    if (strcmp(mode, "error") == 0) {
        printf("%.16e\n", error);
    }
    freeMatrix(mat, ns, ns);
}

/*
 * The following functions are implementations of Forward Euler, Backward
 * Euler, and Crank-Nicolson in one dimension. They all have the same signature
 * so that they can be used interchangeably.
 *
 * Inputs:  up  - pointer to u, the current temperature distribution
 *          r   - ratio involving dx, dt, and diffusivity constant, effectively
 *                  controls how quickly the temperature changes
 *          nx  - # of points in u (length of the array)
 *          nsteps  - # of iterations to run
 *          x0, xL  - boundary conditions
 */

/* Function: FE_1D_Advance
 *
 *  FE_1D_Advance is the naive implementation which loops over each of the
 *  interior points of the rod and calculates the laplacian at each point. The
 *  result is then scaled appropriately by r and added back to u.
 */
void FE_1D_Advance(double **up, int nx, double r, int nsteps, double x0, double xL)
{
    double laplacian[nx];
    laplacian[0] = 0; laplacian[nx-1] = 0;
    for (int i = 0; i < nsteps; i++) {
        for (int j = 1; j < nx-1; j++) {
            laplacian[j] = r * ((*up)[j-1] + (*up)[j+1] - 2 * (*up)[j]);
        }
        vectorAdd(up, laplacian, nx);
    }
}

/* Function: FE_Matrix_Advance
 *
 *  FE_Matrix_Advance bundles the laplacian operator plus the identity as a
 *  matrix so the two operations in the previous implementation can be treated
 *  as a single matrix multiplication.
 */
void FE_Matrix_Advance(double **up, int nx, double r, int nsteps, double x0, double xL)
{
    double **T = createExplicitLaplacian(nx, r);
    for (int i = 0; i < nsteps; i++) {
        matrixVectorMultiply(T, up, nx, nx);
    }
    freeMatrix(T, nx, nx);
}

/* Function: BE_Matrix_Advance
 *
 *  BE_Matrix_Advance takes a similar approach to FE_Matrix_Advance in that it
 *  poses the problem in terms of matrix algebra, but since Backward Euler is an
 *  implicit scheme, the system of equations must be solved at every step. Also
 *  note that the boundary conditions can be affected when essentially
 *  inverting T, whereas in FE, the structure of the matrix maintains them.
 */
void BE_Matrix_Advance(double **up, int nx, double r, int nsteps, double x0, double xL)
{
    double **T = createImplicitLaplacian(nx, nx);

    for (int i = 0; i < nsteps; i++) {
        solveTriDiag(T, up, nx);
        (*up)[0] = x0;
        (*up)[nx-1] = xL;
    }
    freeMatrix(T, nx, nx);
}

/* Function: CN_Matrix_Advance
 *
 *  CN_Matrix_Advance is essentially a combination of Forward and Backward
 *  Euler. Because it takes into account both the explicit and implicit
 *  laplacian, r must be scaled down.
 */
void CN_Matrix_Advance(double **up, int nx, double r, int nsteps, double x0, double xL)
{
    r = r/2;
    double **S = createExplicitLaplacian(nx, nx);
    double **T = createImplicitLaplacian(nx, nx);

    for (int i = 0; i < nsteps; i++) {
        matrixVectorMultiply(S, up, nx, nx);
        solveTriDiag(T, up, nx);
        (*up)[0] = x0;
        (*up)[nx-1] = xL;
    }
    freeMatrix(S,nx,nx);
    freeMatrix(T,nx,nx);
}

/*
 * The two dimensional versions of the time stepping scheme take essentially
 * the same parameters with the exception of a pointer to a matrix rather than
 * an array and additional boundary conditions.
 */

/* Function: FE_2D_Advance_Simple
 *
 *  This version of Forward Euler uses a five point stencil to calculate the
 *  laplacian at every point.
 */
void FE_2D_Advance_Simple(double ***mp, int ns, double r, int nsteps,
        double x0, double xL, double y0, double yL)
{
    double **laplacian = createZeroMatrix(ns, ns);
    for (int t = 0; t < nsteps; t++) {
        for (int i = 1; i < ns -1; i++) {
            for (int j = 1; j < ns-1; j++) {
                laplacian[i][j] = r * ((*mp)[i-1][j] + (*mp)[i+1][j]
                                + (*mp)[i][j-1] + (*mp)[i][j+1] - 4 * (*mp)[i][j]);
            }
        }
        matrixAdd(mp, laplacian, ns, ns);
    }
    freeMatrix(laplacian, ns, ns);
}

/* Function: FE_2D_Advance
 *
 *  FE_2D_Advance takes the vectorized approach to the problem by first
 *  updating the grid using the second spatial derivative in the x direction
 *  (row wise) and then the y direction (column wise).
 */
void FE_2D_Advance(double ***mp, int ns, double r, int nsteps,
        double x0, double xL, double y0, double yL)
{
    double **T = createExplicitLaplacian(ns, r);
    double **m_T, **m_T_T;

    for (int t = 0; t < nsteps; t++) {
        for (int i = 0; i < ns; i++) {
            matrixVectorMultiply(T, (*mp)+i, ns, ns);
        }
        m_T = transposeMatrix((*mp), ns, ns);
        for (int i = 0; i < ns; i++) {
            matrixVectorMultiply(T, m_T+i, ns, ns);
        }
        m_T_T = transposeMatrix(m_T, ns, ns);
        freeMatrix(m_T, ns, ns);
        freeMatrix(*mp, ns, ns);
        *mp = m_T_T;
    }
    freeMatrix(T, ns, ns);
}

/* Function: BE_2D_Advance
 *
 *  BE_2D_Advance does implicitly what FE_2D_Advance does explicitly, taking
 *  the same vectorized approach to the grid.
 */
void BE_2D_Advance(double ***mp, int ns, double r, int nsteps,
        double x0, double xL, double y0, double yL)
{
    double **T = createImplicitLaplacian(ns, r);
    double **m_T, **m_T_T;

    for (int t = 0; t < nsteps; t++) {
        for (int i = 0; i < ns; i++) {
            solveTriDiag(T, (*mp)+i, ns);
        }
        m_T = transposeMatrix((*mp), ns, ns);
        for (int i = 0; i < ns; i++) {
            solveTriDiag(T, m_T+i, ns);
        }
        m_T_T = transposeMatrix(m_T, ns, ns);
        freeMatrix(m_T, ns, ns);
        freeMatrix(*mp, ns, ns);
        *mp = m_T_T;
        /* Enforce BC's */
        for (int i = 0; i < ns; i++) {
            (*mp)[i][0] = y0;
            (*mp)[i][ns-1] = yL;
            (*mp)[0][i] = x0;
            (*mp)[ns-1][i] = xL;
        }
    }
    freeMatrix(T, ns, ns);
}

/* Function: CN_2D_Advance
 *
 *  CN_2D_Advance essentially performs the combination of FE_2D_Advance
 *  followed by BE_Matrix_Advance. Note that the r is cut in half to
 *  accommodate doing both operatios.
 */
void CN_2D_Advance(double ***mp, int ns, double r, int nsteps,
        double x0, double xL, double y0, double yL)
{
    r = r/2;
    double **T = createImplicitLaplacian(ns, r);
    double **S = createExplicitLaplacian(ns, r);
    double **m_T, **m_T_T;
    for (int t = 0; t < nsteps; t++) {
        for (int i = 0; i < ns; i++) {
            matrixVectorMultiply(S, (*mp)+i, ns, ns);
        }
        m_T = transposeMatrix((*mp), ns, ns);
        for (int i = 0; i < ns; i++) {
            matrixVectorMultiply(S, m_T+i, ns, ns);
        }
        m_T_T = transposeMatrix(m_T, ns, ns);
        freeMatrix(*mp, ns, ns);
        *mp = m_T_T;
        for (int i = 0; i < ns; i++) {
            solveTriDiag(T, (*mp)+i, ns);
        }
        m_T = transposeMatrix((*mp), ns, ns);
        for (int i = 0; i < ns; i++) {
            solveTriDiag(T, m_T+i, ns);
        }
        m_T_T = transposeMatrix(m_T, ns, ns);
        freeMatrix(*mp, ns, ns);
        *mp = m_T_T;
        freeMatrix(m_T, ns, ns);

        /* Enforce BC's */
        for (int i = 0; i < ns; i++) {
            (*mp)[i][0] = y0;
            (*mp)[i][ns-1] = yL;
            (*mp)[0][i] = x0;
            (*mp)[ns-1][i] = xL;
        }
    }
    freeMatrix(T, ns, ns);
    freeMatrix(S, ns, ns);
}

/* Function: CN_2D_ADI_Advance
 *
 *  CN_2D_ADI_Advance switches up the order of operations so instead of applying
 *  Forward Euler twice and then Backward Euler twice, it interweaves the two
 *  such that first the x component of the laplacian is explicit and the y
 *  component is implicit and then vice versa.
 */
void CN_2D_ADI_Advance(double ***mp, int ns, double r, int nsteps,
        double x0, double xL, double y0, double yL)
{
    r = r/2;
    double **T = createImplicitLaplacian(ns, r);
    double **S = createExplicitLaplacian(ns, r);
    double **m_T, **m_T_T;

    for (int t = 0; t < nsteps; t++) {
        for (int i = 0; i < ns; i++) {
            matrixVectorMultiply(S, (*mp)+i, ns, ns);
        }
        m_T = transposeMatrix((*mp), ns, ns);
        for (int i = 0; i < ns; i++) {
            solveTriDiag(T, m_T+i, ns);
        }
        m_T_T = transposeMatrix(m_T, ns, ns);
        freeMatrix(*mp, ns, ns);
        *mp = m_T_T;
        for (int i = 0; i < ns; i++) {
            matrixVectorMultiply(S, (*mp)+i, ns, 1);
        }
        m_T = transposeMatrix((*mp), ns, ns);
        for (int i = 0; i < ns; i++) {
            solveTriDiag(T, m_T+i, ns);
        }
        m_T_T = transposeMatrix(m_T, ns, ns);
        freeMatrix(*mp, ns, ns);
        *mp = m_T_T;
        freeMatrix(m_T, ns, ns);

        /* Enforce BC's */
        for (int i = 0; i < ns; i++) {
            (*mp)[i][0] = y0;
            (*mp)[i][ns-1] = yL;
            (*mp)[0][i] = x0;
            (*mp)[ns-1][i] = xL;
        }
    }
    freeMatrix(T, ns, ns);
    freeMatrix(S, ns, ns);
}

/* Function: L2Norm
 *
 * Inputs:
 *      v, exact - vectors to compare distances between
 *      dim      - length of v and exact
 *
 *  Iterates through both arrays and returns the sum of squares of differences
 *  between the two.
 */
double L2Norm(double *v, double *exact, int dim)
{
    double error = 0;
    double calculated, actual;
    for (int i = 0; i < dim; i++) {
            calculated = v[i];
            actual = exact[i];
            error += simple_square(calculated - actual);
    }
    return error;
}

/* Function: printMatrix
 *
 * Inputs:
 *      A   - grid of doubles
 *      n   - # of rows
 *      m   - # of columns
 *      nd  - # of digits of precision
 *
 *  Given a matrix and its dimensions, prints its contents out to stdout
 *  with the given precision.
 */
void printMatrix(double **A, int n, int m, int nd)
{
    for (int i = 0; i < n; i++){
        printVec(A[i], m, nd);
    }
    printf("\n\n");
}

/* Overloaded function for default precision of 5 */
void printMatrix(double **A, int n, int m)
{
    for (int i = 0; i < n; i++){
        printVec(A[i], m);
    }
    printf("\n\n");
}

/* Function: linspace
 *
 * Inputs:
 *      start   - starting value of the sequence
 *      stop    - ending value of the sequence
 *      n       - # of values in the sequence
 *
 *  Creates a dynamically allocated array of size n populated with a evenly
 *  spaced values between start and stop.
 */
double *linspace(double start, double stop, int n)
{
    double *result = (double *) malloc(n * sizeof(double));
    double dx = (stop - start)/(n-1);
    for (int i = 0; i < n; i++) {
        result[i] = start + i * dx;
    }
    return result;
}

/* Function: vecApply
 *
 * Inputs:
 *      f   - generic function to call for each of the elements
 *      arr - array of doubles (to be passed into f)
 *      len - length of the array
 *
 *  vecApply returns a dynamically allocated copy of f applied to each of the
 *  elements in arr. The contents of arr are not changed.
 */
double *vecApply(double (*f)(double), double *arr, int len)
{
    double *result = (double *) malloc(len * sizeof(double));
    for (int i = 0; i < len; i++) {
        result[i] = f(arr[i]);
    }
    return result;
}

/* Function: transposeMatrix
 *
 * Inputs:
 *      A   - the original matrix
 *      n,m - # of rows and columns
 *
 *  Returns a dynamically allocated copy of the transpose of A. The contents of
 *  A are not modified.
 */
double **transposeMatrix(double **A, int n, int m)
{
    double **A_T = createZeroMatrix(m, n);
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < n; col++) {
            A_T[row][col] = A[col][row];
        }
    }
    return A_T;
}

/* Function: matrixMultiply
 *
 * Inputs:
 *      A, B  - matrices to multiply
 *      size  - dimension of A,B
 *
 *  Assumes that A and B are square matrices and returns a dymanically allocated
 *  copy of their product.
 */
double **matrixMultiply(double **A, double **B, int size)
{
    /* First transpose the matrix so we can express the multiplication as a dot
     * product, saving cache misses. */
    double **B_T = transposeMatrix(B, size, size);
    double **result = createZeroMatrix(size, size);

    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            for (int i = 0; i < size; i++){
                result[row][col] += A[row][i] * B_T[col][i];
            }
        }
    }
    freeMatrix(B_T, size, size);
    return result;
}

/* Function: matrixAdd
 *
 * Inputs:
 *      Ap  - pointer to matrix representing A
 *      B   - matrix representing B
 *      n,m - dimensions of A and B
 *
 *  Computes the matrix sum of A and B. Note that A is modified in place!
 */
void matrixAdd(double ***Ap, double **B, int n, int m)
{
    for (int i = 0; i < n; i++) {
        vectorAdd((*Ap) + i, B[i], m);
    }
}

/* Function: matrixVectorMultiply
 *
 * Inputs:
 *      T   - matrix representing T
 *      vp  - pointer to a vector v
 *      n,m - dimensions of T
 *
 *  Computes the product of T applied to v, which is modified in place.
 */
void matrixVectorMultiply(double **T, double **vp, int n, int m)
{
    double temp[m];
    memset(temp, 0, m * sizeof(double));
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < m; col++) {
            temp[row] += T[row][col] * (*vp)[col];
        }
    }
    for (int i = 0; i < m; i++) {
        (*vp)[i] = temp[i];
    }
}

/* Function: solveTriDiag
 *
 * Inputs:
 *      T   - matrix representing T
 *      up  - pointer to a vector u
 *      n   - length of u, # of rows/columns in T
 *
 *  Applies Thomas' algorithm for solving tridiagonal matrices. The primary
 *  use of this function is for solving the system of equations that arises
 *  from Backward Euler.
 */
void solveTriDiag(double **T, double **up, int n)
{
    double a[n], b[n], c[n], c_prime[n], d[n], u_new[n];
    for (int i = 0; i < n; i++) {
        b[i] = T[i][i];
    }
    for (int i = 0; i < n-1; i++) {
        a[i+1] = T[i+1][i];
        c[i] = T[i][i+1];
    }
    c_prime[0] = c[0]/b[0];
    u_new[0] = (*up)[0]/b[0];
    for (int i = 1; i <  n; i++) {
        c_prime[i] = c[i] / (b[i]  - c_prime[i-1] * a[i]);
        u_new[i] = ((*up)[i] - u_new[i-1] * a[i])/(b[i] - c_prime[i-1] * a[i]);
    }

    for (int i = n-2; i >= 0; i--) {
        u_new[i] -= c_prime[i] * u_new[i+1];
    }

    for (int i = 0; i < n; i++) {
        (*up)[i] = u_new[i];
    }
}

/* Helper function which returns an n by m matrix filled with 0's */
double **createZeroMatrix(int n, int m)
{
    double **result = (double **) malloc(n * sizeof(double*));
    for (int row = 0; row < n; row++)
        result[row] = (double *) calloc(m , sizeof(double));
    return result;
}

/* Returns the matrix representation of 1D Forward Euler. */
double **createExplicitLaplacian(int dim, double r)
{
    double **result = createZeroMatrix(dim, dim);
    for (int i = 0; i < dim; i++) {
        if (i == 0 || i == dim - 1) {
            result[i][i] = 1;
        } else {
            result[i][i] = 1 - 2 * r;
            result[i][i-1] = r;
            result[i][i+1] = r;
        }
    }
    return result;
}
/* Returns the matrix representation of 1D Backward Euler. */
double **createImplicitLaplacian(int dim, double r)
{
    double **result = createZeroMatrix(dim, dim);
    for (int i = 0; i < dim; i++) {
        if (i == 0 || i == dim - 1) {
            result[i][i] = 1;
        } else {
            result[i][i] = 1 + 2 * r;
            result[i][i-1] = -r;
            result[i][i+1] = -r;
        }
    }
    return result;
}

/* Handles memory deallocation for an n x m matrix */
void freeMatrix(double **arr, int n, int m)
{
    for (int row = 0; row < n; row++)
        free(arr[row]);
    free(arr);
}

/* Scales every element in a vector u by c in place. */
void vectorScale(double **up, double c, int len)
{
    for (int i = 0; i < len; i++) {
        (*up)[i] *= c;
    }
}

/* Computes the vector sum of u and v and stores the result in u. */
void vectorAdd(double **up, double *v, int len)
{
    for (int i = 0; i < len; i++) {
        (*up)[i] += v[i];
    }
}

/* Prints the contents of an array v with nd digits of precision. */
void printVec(double *v, int n, int nd)
{
    for (int i = 0; i < n-1; i++){
        printf("%*.*g ", nd+2, nd, v[i]);
    }
    printf("%*.*g\n", nd+2, nd, v[n-1]);
}
/* Overloaded method for a default nd of 5. */
void printVec(double *v, int n)
{
    printVec(v, n, 5);
}

double simple_square(double x) { return x * x; }

/* Used for debugging purposes */
void validateRoutines()
{
    int N = 11;
    printf("Testing linspace function\n");
    double *x_vals = linspace(0, 10, N);
    printf("xvals: ");
    printVec(x_vals, N);
    printf("Testing vecApply\n");
    double *y_vals = vecApply(simple_square, x_vals, N);
    printf("yvals: ");
    printVec(y_vals, N);

    N = 3;
    double **A = createZeroMatrix(N,N);
    double **B = createZeroMatrix(N,N);
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            A[i][j] = 3 * i + j;
            B[i][j] = 3 * i + j;
        }
    }
    printf("Used createZeroMatrix to generate skeletons for A and B\n\n");

    printf("Testing matrix multiplcation\n");
    double **res = matrixMultiply(A, B, N);
    for (int i = 0; i < 3; i++){
        printf("%f %f %f\n", res[i][0],res[i][1],res[i][2]);
    }
    printf("\n\n");


    double *v = (double *) malloc(3 * sizeof(double));
    printVec(v, 3);
    v[0] = 1;
    v[1] = 2;
    v[2] = 3;
    printVec(v, 3);
    printf("A:\n");
    for (int i = 0; i < 3; i++){
        printf("%f %f %f\n", A[i][0],A[i][1],A[i][2]);
    }
    printf("v:\n");
    printVec(v, 3);
    printf("\n\n");
    matrixVectorMultiply(A, &v, 3, 3);
    printf("Testing matrix vector multiplication..\n");
    printVec(v, 3);

    printf("Testing matrix addition\n");
    matrixAdd(&A, B, 3, 3);
    printf("A:\n");
    for (int i = 0; i < 3; i++){
        printf("%f %f %f\n", A[i][0],A[i][1],A[i][2]);
    }
    printf("\n\n");


    /* Cleanup */
    free(y_vals);
    free(v);
    freeMatrix(A, 3, 3);
    freeMatrix(B, 3, 3);
}
