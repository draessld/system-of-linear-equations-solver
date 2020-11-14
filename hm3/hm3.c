#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <argp.h>
#include <stdbool.h>
#include <sys/time.h>

#define EPSILON 1e-5

size_t n;   //  matrix/vector size
int k = 10; //  max steps
size_t z;   //  cr number values

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  parser part
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//  parsing setting
const char *argp_program_version = "1.0";
const char *argp_program_bug_address = "<draesdom@cvut.cz>";
static char doc[] = "Solve algebraic equations";
static char args_doc[] = "[FILENAMES]...";
static struct argp_option options[] = {
    {"2d-array", 'a', 0, 0, "Save matrix as 2d array format."},
    {"compess-row", 'r', 0, 0, "Save matrix as compress row."},
    {"gradient-descentd", 'g', 0, 0, "Use gradient descent method."},
    {"conjugate-gradient", 'c', 0, 0, "Use conjugate gradient method."},
    {"output", 'o', "OUTFILE", 0, "Output file path."},
    {"size", 'n', "N", 0, "squared matrix size"},
    {"number-of-steps", 'k', "K", 0, "maximal number of iterative steps"},
    {0}};

struct arguments
{
    enum method
    {
        METHOD_GD,
        METHOD_CG,
        NO_METHOD
    } method;
    enum format
    {
        FORMAT_2D,
        FORMAT_CR,
        NO_FORMAT
    } format;
    bool count_n;
    char *outfile;
};

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
        //  check if format already setted
        if (arguments->format != NO_FORMAT)
        {
            char c = 'c';
            printf("matrix format already setted on %s. Replace it (y/n):", format_to_string(arguments->format));
            while (((c != 'n') && (c != 'y')))
            {
                scanf(" %c", &c);
            }
            if (c == 'y')
            {
                arguments->method = FORMAT_2D;
                printf("replacing on 2D array format.\n");
            }
        }
        else
            arguments->format = FORMAT_2D;
        break;
    case 'r':
        //  check if format already setted
        if (arguments->format != NO_FORMAT)
        {
            char c = 'c';
            printf("matrix format already setted on %s. Replace it (y/n):", format_to_string(arguments->format));
            while (((c != 'n') && (c != 'y')))
            {
                scanf(" %c", &c);
            }
            if (c == 'y')
            {
                arguments->method = FORMAT_CR;
                printf("replacing on compress row format.\n");
            }
        }
        else
            arguments->format = FORMAT_CR;
        break;
    //  method keys
    case 'g':
        //  check if method already setted
        if (arguments->method != NO_METHOD)
        {
            char c;
            printf("method already setted on %s. Replace it (y/n):", method_to_string(arguments->method));
            while (((c != 'n') && (c != 'y')))
            {
                scanf(" %c", &c);
            }
            if (c == 'y')
            {
                arguments->method = METHOD_GD;
                printf("replacing on gradient descent method.\n");
            }
        }
        else
            arguments->method = METHOD_GD;
        break;
    case 'c':
        //  check if method already setted
        if (arguments->method != NO_METHOD)
        {
            char c;
            printf("method already setted on %s. Replace it (y/n):", method_to_string(arguments->method));
            while (((c != 'n') && (c != 'y')))
            {
                scanf(" %c", &c);
            }
            if (c == 'y')
            {
                arguments->method = METHOD_CG;
                printf("replacing on gradient descent method.\n");
            }
        }
        else
        {
            arguments->method = METHOD_CG;
        }
        break;
    //  output files
    case 'o':
        arguments->outfile = arg;
        break;
    //  other
    case 'n':
        arguments->count_n = false;
        n = atoi(arg);
        break;
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

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  printers
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//  control print of matrix
void print_m(double **a)
{

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            printf("%lf ", a[i][j]);
        printf("\n");
    }
}

//  control print of vector
void print_v(double *b)
{
    for (int j = 0; j < n; j++)
        printf("%lf ", b[j]);
    printf("\n");
}

//  control print csr matrix
void print_csr(double **a, int *adr)
{
    printf("values: ");
    for (int i = 0; i < z; i++)
    {
        printf("%-5lf ", a[1][i]);
    }
    printf("\ncol: ");

    for (int i = 0; i < z; i++)
    {
        printf("%-5lf ", a[0][i]);
    }

    printf("\nrow_adr: ");
    for (int i = 0; i < n; ++i)
        printf("%-5d ", adr[i]);

    printf("\n");
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  file read/write functions
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//  read matrix as 2d
int read_2d_matrix(double **a, const char *filename)
{

    FILE *pf;
    pf = fopen(filename, "r");
    if (pf == NULL)
        return 0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            fscanf(pf, "%lf", a[i] + j);
    }

    fclose(pf);
    return 1;
}

//  read matrix as compressed row
int read_csr_matrix(double **val, int *adr, const char *filename)
{
    FILE *pf;
    pf = fopen(filename, "r");
    if (pf == NULL)
        return 0;

    fscanf(pf, "%ld %ld", &n, &z);
    // printf("%ld, %ld\n", n, z);
    val[0] = (double *)malloc(z * sizeof(double)); //  column index
    val[1] = (double *)malloc(z * sizeof(double)); //  value
    adr = (int *)realloc(adr, n);                  //  index of first row value

    //  first read
    int d = 0;
    int old_d = 0;
    int adr_i = 0;

    //  next read
    for (int i = 0; i < z; ++i)
    {
        fscanf(pf, " %d %lf %lf", &d, &val[0][i], &val[1][i]);
        // printf("%d na indexu %d - %d\n", d,i, d - old_d);
        if (d - old_d == 1)
        { //  new value on next line
            adr[adr_i] = i;
            // printf(" pisu %d na index %d\n", i,adr_i);
            old_d = d;
            adr_i++;
        }
        else if (d - old_d > 1)
        { //    d-old_d zero rows - write -1 and move
            for (int j = 0; j < d - old_d - 1; j++)
            {
                adr[adr_i] = -1;
                // printf(" pisu -1 na index %d\n", adr_i);
                adr_i++;
            }
            // printf(" pisu %d na index %d\n", i,adr_i);
            adr[adr_i] = i;
            adr_i++;
            old_d = d;
        }
    }

    fclose(pf);
    return 1;
}

//  read vector
int read_vector(double *a, const char *filename)
{
    FILE *pf;
    pf = fopen(filename, "r");
    if (pf == NULL)
        return 0;

    for (int i = 0; i < n; ++i)
        fscanf(pf, "%lf\n", &a[i]);

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
        fprintf(pf, "%lf ", a[i]);

    fclose(pf);
    return 1;
}

//  get number of equations by size of vector
int read_n(const char *filename)
{
    double x;

    FILE *pf;
    pf = fopen(filename, "r");
    if (pf == NULL)
        return 0;

    while (fscanf(pf, " %lf", &x) != EOF)
    {
        n++;
    }

    fclose(pf);
    return 1;
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  basic calculations
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//  matrix * column vector
void multiply_m_v(double **a, double *b, double *res)
{
    for (int i = 0; i < n; i++)
    {
        res[i] = 0;
        for (int j = 0; j < n; j++)
        {
            // printf("%lf + (%lf * %lf) = %lf + %lf  ",res[j],a[j][i], b[i],res[j], a[j][i] * b[i]);
            res[i] += (a[i][j] * b[j]);
            // printf(" = %lf; \n",res[j]);
        }
    }
}

void csr_multiply_m_v(double **a, int *rowPtr, double *b, double *res)
{
    for (int i = 0; i < n; i++)
    {
        res[i] = 0;
        if (rowPtr[i] == -1)
            res[i] = 0;
        else
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; k++)
                res[i] += a[1][k] * b[(int)a[0][k]];
    }
}

//  substract two vectors
void subtract_v_v(double *a, double *b, double *res)
{
    for (int i = 0; i < n; i++)
        res[i] = a[i] - b[i];
}

//  vector zvetseny o konstantu
void multiply_v_c(double *a, double c, double *res)
{
    for (int i = 0; i < n; i++)
        res[i] = a[i] * c;
}

//  create copy of vector
void copy_vector(double *s, double *r)
{
    for (int i = 0; i < n; i++)
        r[i] = s[i];
}

//  additive two vector
void addition_v_v(double *a, double *b, double *res)
{
    for (int i = 0; i < n; i++)
        res[i] = a[i] + b[i];
}

//  sloupcovy vector * radkovy vector
double multiply_v_v(double *a, double *b)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += a[i] * b[i];
    return sum;
}

//  check if vector = zero vector
bool empty(double *x)
{
    for (int i = 0; i < n; i++)
        if (x[i] != 0)
            return 1;
    return 0;
}

//  return vector size
double vector_size(double *vector)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += pow(vector[i], 2.);
    return sqrt(sum);
}

//  return size
double gen_vector_size(double *a, double **A)
{
    double r[n];
    multiply_m_v(A, a, r);
    return multiply_v_v(a, r);
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  methods
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//  gradient descent method
void gradient_descent_method(double **A, double *b, double *xk)
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
    multiply_m_v(A, xk, Ax);
    subtract_v_v(b, Ax, rk);

    double rho = multiply_v_v(rk, rk);
    if (sqrt(rho) == 0)
        return;

    //  lets do it! :))
    for (int i = 0; i < k; i++)
    {
        //  ak = rk'*rk / rk'*A*rk
        multiply_m_v(A, rk, Ax);
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
void conjugate_gradient_method(double **A, double *b, double *xk)
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
    multiply_m_v(A, xk, rk);
    subtract_v_v(b, rk, rk);

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
        multiply_m_v(A, sk, pom1);
        ak = rho / multiply_v_v(sk, pom1);

        //  rk_1 = rk - ak*A*sk
        multiply_v_c(pom1, ak, pom1);
        subtract_v_v(rk, pom1, rk);

        //  xk_1 = xk + ak*sk
        multiply_v_c(sk, ak, pom1);
        addition_v_v(xk, pom1, xk);

        //  time to end?
        if (gen_vector_size(rk, A) < EPSILON)
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

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  main
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

int main(int argc, char *argv[])
{
    struct arguments arguments;

    //  default arguments
    arguments.format = NO_FORMAT; //  default 2d array atric format
    arguments.method = NO_METHOD; //  default gradient descent method
    arguments.count_n = true;     //  default unknown size of matrix - calculate it while reading
    arguments.outfile = NULL;     //  default output file - no file

    //  input filenames
    char *matrix_input = argv[argc - 2];
    char *vector_input = argv[argc - 1];

    //  get arguments
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    //  read matrix size if required
    if (arguments.count_n)
        if (read_n(vector_input) == 0)
        {
            printf("filename %s does not exist", vector_input);
            return 1;
        }

    //  print basic information about calculation
    printf("There are %ld algebraic equations to solve\n", n);
    printf("Procces terminated after %d steps\n", k);
    printf("Matrix format: %s.\n", format_to_string(arguments.format));
    printf("Method: %s.\n", method_to_string(arguments.method));
    printf("INPUTS READING:\n");
    printf("------------------------------------------\n");

    //  declare matrix
    double **matrix;
    int *adr;

    // declare vector and allocate mem space
    double *vector = (double *)malloc(n * sizeof(double));

    //  read matrix
    switch (arguments.format)
    {
    case FORMAT_2D:
        //  allocate mem space for 2d array
        matrix = (double **)malloc(n * sizeof(double *));
        for (int i = 0; i < n; i++)
            matrix[i] = (double *)malloc(n * sizeof(double));
        if (read_2d_matrix(matrix, matrix_input) == 0)
        {
            printf("filename %s does not exist", matrix_input);
            return 1;
        }
        break;
    case FORMAT_CR:
        //  allocate mem space for 2d array
        matrix = (double **)malloc(2 * sizeof(double *));
        adr = (int *)malloc(1 * sizeof(int));
        if (read_csr_matrix(matrix, adr, matrix_input) == 0)
        {
            printf("filename %s does not exist", matrix_input);
            return 1;
        }
        break;
    default:
        printf("Unknown format - Exiting\n");
        return 1;
        break;
    }

    //  log inform
    printf("Matrix successfuly readed.\n");

    // read vector
    if (read_vector(vector, vector_input) == 0)
    {
        printf("filename %s does not exist", vector_input);
        return 1;
    }

    //  log inform
    printf("Vector successfuly readed.\n");

    // print_csr(matrix, adr);
    // print_v(vector);

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
    // switch (arguments.method)
    // {
    // case METHOD_GD:
    //     // gradient_descent_method(matrix, vector, output);
    //     break;
    // case METHOD_CG:
    //     // conjugate_gradient_method(matrix, vector, output);
    // default:
    //     break;
    // }

    csr_multiply_m_v(matrix, adr, vector, output);

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

    gettimeofday(&stop, NULL);
    printf("time %lu us\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);

    //  tidy up
    if (arguments.format == FORMAT_2D)
    {
        for (int i = 0; i < n; i++)
            free(matrix[i]);
        free(matrix);
    }
    free(vector);
    // free(output);

    return 0;
}
