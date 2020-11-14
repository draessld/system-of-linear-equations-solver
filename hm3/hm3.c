#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <argp.h>
#include <stdbool.h>
#include <sys/time.h>
#include <assert.h>

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *  globals and constants
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

#define EPSILON 1e-5

size_t n;    //  matrix/vector size
size_t z;    //  cr number values
int k = 100; //  max steps

double **matrix;
double *vector;
int *rowPtr;

typedef struct arguments
{
    enum method
    {
        METHOD_GD,
        METHOD_CG
    } method;
    enum format
    {
        FORMAT_2D,
        FORMAT_CR
    } format;
    char *outfile;
}Arguments;

Arguments arguments;

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *  parser setting
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

//  parsing setting
const char *argp_program_version = "1.0";
const char *argp_program_bug_address = "<draesdom@cvut.cz>";
static char doc[] = "Solve algebraic equations";
static char args_doc[] = "[FILENAMES]...";
static struct argp_option options[] = {
    {"2d-array", 'a', 0, 0, "Save matrix as 2d array format."},
    {"compess-row", 'r', 0, 0, "Save matrix as compress row."},
    {"gradient-descent", 'g', 0, 0, "Use gradient descent method."},
    {"conjugate-gradient", 'c', 0, 0, "Use conjugate gradient method."},
    {"output", 'o', "OUTFILE", 0, "Output file path."},
    {"number-of-steps", 'k', "K", 0, "maximal number of iterative steps"},
    {0}};

//  return string name of enum format
char *format_to_string(enum format a)
{
    switch (a)
    {
    case FORMAT_2D:
        return "2d_array";
        break;
    case FORMAT_CR:
        return "compress row";
        break;
    default:
        return "unknown_format";
        break;
    }
}

//  return string name of enum method
char *method_to_string(enum method a)
{
    switch (a)
    {
    case METHOD_GD:
        return "gradient descent method";
        break;
    case METHOD_CG:
        return "conjugate gradient method";
        break;
    default:
        return "unknown_method";
        break;
    }
}

//  set method, format and other program options
static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;
    switch (key)
    {
    //  matrix format keys
    case 'a':
        arguments->format = FORMAT_2D;
        break;
    case 'r':
        arguments->format = FORMAT_CR;
        break;
    //  method keys
    case 'g':
        arguments->method = METHOD_GD;
        break;
    case 'c':
        arguments->method = METHOD_CG;
        break;
    //  output files
    case 'o':
        arguments->outfile = arg;
        break;
    //  other
    case 'k':
        k = atoi(arg);
        break;
    case ARGP_KEY_ARG:
        return 0;
    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *  readers/writers
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

//  read matrix
int read_matrix(const char *filename)
{
    FILE *pf;
    pf = fopen(filename, "r");
    if (pf == NULL)
        return 0;

    switch (arguments.format)
    {
    case FORMAT_2D:
        printf("Reading 2d array format\n");
        fscanf(pf, " %ld", &n);
        //  read like 2d array
        matrix = calloc(n,sizeof(double*)); //  column index
        for (int i = 0; i < n; ++i)
        {
            matrix[i] = calloc(n,sizeof(double));
            for (int j = 0; j < n; ++j)
                fscanf(pf, "%lf", matrix[i] + j);
        }
        break;
    case FORMAT_CR:
        printf("Reading CSR format\n");
        fscanf(pf, "%ld %ld", &n, &z);
        matrix = calloc(2,sizeof(double*));
        matrix[0] = (double *)malloc(z * sizeof(double)); //  column index
        matrix[1] = (double *)malloc(z * sizeof(double)); //  values
        rowPtr = (int *)malloc((n + 1) * sizeof(int));  //  index of first row value

        int d = 0;
        int old_d = 0;
        int rowPtr_i = 0;

        //  next read
        for (int i = 0; i < z; ++i)
        {
            fscanf(pf, " %d %lf %lf", &d, &matrix[0][i], &matrix[1][i]);
            if (d - old_d == 1)
            { //  new value on next line
                rowPtr[rowPtr_i] = i;
                old_d = d;
                rowPtr_i++;
            }
            else if (d - old_d > 1)
            { //    d-old_d zero rows - write -1 and move
                for (int j = 0; j < d - old_d - 1; j++)
                {
                    rowPtr[rowPtr_i] = -1;
                    rowPtr_i++;
                }
                rowPtr[rowPtr_i] = i;
                rowPtr_i++;
                old_d = d;
            }
        }
        rowPtr[n] = z;
        break;
    default:
        printf("Unknown format!\n");
        return -1;
        break;
    }

    fclose(pf);
    return 1;
}

//  read vector
int read_vector(const char *filename)
{
    FILE *pf;
    pf = fopen(filename, "r");
    if (pf == NULL)
        return 0;

    for (int i = 0; i < n; ++i)
        fscanf(pf, "%lf\n", &vector[i]);

    fclose(pf);
    return 1;
}

//  write output
int write_output(double *a, const char *filename)
{
    FILE *pf;
    pf = fopen(filename, "w");
    if (pf == NULL)
        return 0;

    for (int i = 0; i < n; ++i)
        fprintf(pf, "%1.15e\n", a[i]);

    fclose(pf);
    return 1;
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *  calculations
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

//  matric * vector
void multiply_m_v(double *v, double *res)
{
    switch (arguments.format)
    {
    case FORMAT_2D:
        for (int i = 0; i < n; i++)
        {
            res[i] = 0;
            for (int j = 0; j < n; j++)
                res[i] += (matrix[i][j] * v[j]);
        }
        break;
    case FORMAT_CR:
        for (int i = 0; i < n; i++)
        {
            res[i] = 0;
            if (rowPtr[i] == -1)
                res[i] = 0;
            else
                for (int k = rowPtr[i]; k < rowPtr[i + 1]; k++)
                    res[i] += matrix[1][k] * v[(int)matrix[0][k] - 1];
        }
        break;
    default:
        printf("Unknown format\n");
        break;
    }
}

//  vector * vector
double multiply_v_v(double *a, double *b)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += a[i] * b[i];
    return sum;
}

//  vector * constant
void multiply_v_c(double *a, double c, double *res)
{
    for (int i = 0; i < n; i++)
        res[i] = a[i] * c;
}

//  vector - vector
void subtract_v_v(double *a, double *b, double *res)
{
    for (int i = 0; i < n; i++)
        res[i] = a[i] - b[i];
}

//  vector + vector
void addition_v_v(double *a, double *b, double *res)
{
    for (int i = 0; i < n; i++)
        res[i] = a[i] + b[i];
}

//  return ordinary vector size
double vector_size(double *v)
{
    return sqrt(multiply_v_v(v, v));
}

//  return vector size in generalized sense
double gen_vector_size(double *v)
{
    double r[n];
    multiply_m_v(v, r);
    return sqrt(multiply_v_v(v, r));
}

//  create copy of vector
void copy_vector(double *s, double *r)
{
    for (int i = 0; i < n; i++)
        r[i] = s[i];
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *  methods
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

//  gradient descent method
void gradient_descent_method(double *xk)
{

    //  initialization vector x
    double Ax[n];
    double rk[n];
    double ar[n];
    for (int i = 0; i < n; i++)
    {
        xk[i] = 0;
        rk[i] = 0;
        Ax[i] = 0;
        ar[i] = 0;
    }

    double ak = 0; //  lenght of step

    //  rk = b - A*xk
    multiply_m_v(xk, Ax);
    subtract_v_v(vector, Ax, rk);

    double rho = multiply_v_v(rk, rk);
    if (sqrt(rho) == 0)
        return;

    //  lets do it! :))
    for (int i = 0; i < k; i++)
    {
        //  ak = rk'*rk / rk'*A*rk
        multiply_m_v(rk, Ax);
        ak = rho / multiply_v_v(rk, Ax);

        //  xk_1 = xk + ak*rk
        multiply_v_c(rk, ak, ar);
        addition_v_v(xk, ar, xk);

        //  rk = b - A*xk
        multiply_v_c(Ax, ak, Ax);
        subtract_v_v(rk, Ax, rk);

        rho = multiply_v_v(rk, rk);

        if (sqrt(rho) < EPSILON)
            break;
    }
}

//  conjugate gradient method
void conjugate_gradient_method(double *xk)
{
    //  initialization vectors
    double sk[n];
    double pom1[n];
    double rk[n];
    for (int i = 0; i < n; i++)
    {
        xk[i] = 0;
        sk[i] = 0;
        pom1[i] = 0;
        rk[i] = 0;
    }

    // rk = b - A*xk
    multiply_m_v(xk, rk);
    subtract_v_v(vector, rk, rk);

    // sk = rk
    copy_vector(rk, sk);

    double ak = 0; //  lenght of step
    double bk = 0; //  pom constant
    double rho = multiply_v_v(rk, rk);
    if (sqrt(rho) == 0)
        return;

    for (int i = 0; i < k; i++)
    {
        //  ak = rk'*rk / sk'*A*sk
        multiply_m_v(sk, pom1);
        ak = rho / multiply_v_v(sk, pom1);

        //  rk_1 = rk - ak*A*sk
        multiply_v_c(pom1, ak, pom1);
        subtract_v_v(rk, pom1, rk);

        //  xk_1 = xk + ak*sk
        multiply_v_c(sk, ak, pom1);
        addition_v_v(xk, pom1, xk);

        //  time to end?
        if (gen_vector_size(rk) < EPSILON)
            break;

        //  bk = rk_1' * rk_1 / rk' * rk
        bk = rho;
        rho = multiply_v_v(rk, rk);
        bk = rho / bk;

        //  sk_1 = rk_1 + bk*sk
        multiply_v_c(sk, bk, sk);
        addition_v_v(rk, sk, sk);
    }
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *  other
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

//  control print of vector
void print_v(double *b)
{
    for (int j = 0; j < n; j++)
        printf("%1.15e\n", b[j]);
    printf("\n");
}

//  free allocated memory
void tidy_up(double *output){
    if (arguments.format == FORMAT_2D)
    {
        for (int i = 0; i < n; i++)
            free(matrix[i]);
    }
    else
    {
        free(matrix[0]);
        free(matrix[1]);
    free(rowPtr);
    }
    free(matrix);
    free(vector);
    free(output);
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *  main
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

int main(int argc, char *argv[])
{

    //  default arguments
    arguments.outfile = NULL;     //  default output file - no file

    //  input filenames
    char *matrix_input = argv[argc - 2];
    char *vector_input = argv[argc - 1];

    //  get arguments
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    //  print basic information about calculation
    printf("Procces terminated after %d steps\n", k);
    printf("Matrix format: %s.\n", format_to_string(arguments.format));
    printf("Method: %s.\n", method_to_string(arguments.method));
    printf("------------------------------------------\n");
    printf("READING:\n");
    printf("------------------------------------------\n");

    //  read matrix
    assert(read_matrix(matrix_input) == 1);
    printf("Matrix successfuly readed.\n");

    // read vector
    vector = (double *)malloc(n * sizeof(double));
    assert(read_vector(vector_input) == 1);
    printf("Vector successfuly readed.\n");

    //  declare output vector
    double *output = (double *)malloc(n * sizeof(double));

    //  log inform
    printf("------------------------------------------\n");
    printf("CALCULATION:\n");
    printf("------------------------------------------\n");

    //  start time calculation
    struct timeval stop, start;
    gettimeofday(&start, NULL);

    //  switch by method
    switch (arguments.method)
    {
    case METHOD_GD:
        gradient_descent_method(output);
        break;
    case METHOD_CG:
        conjugate_gradient_method(output);
    default:
        break;
    }

    //  log inform
    printf("Calculation done.\n");

    //  print result
    if (arguments.outfile)
    {
        write_output(output, arguments.outfile);
        printf("Result at %s.\n", arguments.outfile);
    }
    else
    {
        printf("Result:\n");
        print_v(output);
    }

    //  stop time
    gettimeofday(&stop, NULL);
    printf("time %lu us\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);

    //  tidy up
    tidy_up(output);

    return 0;
}
