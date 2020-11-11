#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <argp.h>
#include <stdbool.h>
#include <sys/time.h>

#define EPSILON 1e-5

size_t n;   //  matrix/vector size
int k = 10; //  max steps

//  parsing setting
const char *argp_program_version = "1.0";
const char *argp_program_bug_address = "<draesdom@cvut.cz>";
static char doc[] = "Solve algebraic equations";
static char args_doc[] = "[FILENAMES]...";
static struct argp_option options[] = {
    {"2d-array", 'a', 0, 0, "Save matrix as 2d array format."},
    {"gradient-descentd", 'g', 0, 0, "Use gradient descent method."},
    {"conjugate-gradient", 'c', 0, 0, "Use conjugate gradient method."},
    {"output", 'o', "OUTFILE", 0, "Output file path."},
    {"size", 'n', "N", 0, "squared matrix size"},
    {"number--of--steps", 'k', "K", 0, "maximal number of iterative steps"},
    {0}};

struct arguments
{
    enum method
    {
        METHOD_GD,
        METHOD_CG
    } method;
    enum format
    {
        FORMAT_2D
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
        if (arguments->format != FORMAT_2D)
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
        break;
    //  method keys
    case 'g':
        //  check if method already setted
        if (arguments->method != METHOD_GD)
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
        break;
    case 'c':
        //  check if method already setted
        if (arguments->method != METHOD_CG)
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

//  return lenght of vector
double vector_lenght(double *vector)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += pow(vector[i], 2.);
    return sqrt(sum);
}

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

//  read matrix as 2d
int read_matrix(double **a, const char *filename)
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

//  read vector
int read_vector(double *a, const char *filename)
{
    FILE *pf;
    pf = fopen(filename, "r");
    if (pf == NULL)
        return 0;

    for (int i = 0; i < n; ++i)
        fscanf(pf, "%lf", &a[i]);

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

//  matrix * column vector
void muptiply_m_v(double **a, double *b, double *res)
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

//  check if vector is zero
bool empty(double *x)
{
    for (int i = 0; i < n; i++)
        if (x[i] != 0)
            return 1;
    return 0;
}

//  gradient descent method
void gradient_descent_method(double **matrix, double *b, double *x)
{

    //  initialization vector x
    double pom[n];
    double pom1[n];
    for (int i = 0; i < n; i++)
    {
        x[i] = 0;
        pom[i] = 0;
        pom1[i] = 0;
    }

    double *r = (double *)malloc(n * sizeof(double));
    muptiply_m_v(matrix, x, pom);
    subtract_v_v(b, pom, pom);
    if (empty(pom))
        return;

    copy_vector(pom, r);
    double lb = vector_lenght(b);
    double lr = 0;
    double a = 0; //  lenght of step
    double e = 0; //  difference
    bool p = true;

    //  lets do it! :))
    do
    {
        muptiply_m_v(matrix, r, pom);
        if (p)
        {
            p = false;
            a = multiply_v_v(r, r) / multiply_v_v(r, pom);
        }
        else
            a = pow(lr, 2) / multiply_v_v(r, pom);
        multiply_v_c(r, a, pom1);
        addition_v_v(x, pom1, x);
        multiply_v_c(pom, a, pom1);
        subtract_v_v(r, pom1, r);
        lr = vector_lenght(r);
        if (lb == 0)
            e = lr;
        else
            e = lr / lb;
        // printf("e: %lf \n",e);
    } while (e > EPSILON);

    free(r);
}

double gen_vector_lenght(double *a, double **A)
{
    double r[n];
    muptiply_m_v(A, a, r);
    return multiply_v_v(a, r);
}

//  conjugate gradient method
void conjugate_gradient_method(double **A, double *b, double *xk)
{
    //  initialization vectors
    double sk[n];
    double pom1[n];
    double rk[n];
    // double *rk = (double *)malloc(n * sizeof(double));
    // double *rk_1 = (double *)malloc(n * sizeof(double));
    
    for (int i = 0; i < n; i++)
    {
        xk[i] = 0;
        sk[i] = 0;
        pom1[i] = 0;
        rk[i] = 0;
        // rk_1[i] = 0;
    }

    // rk = b - A*xk
    muptiply_m_v(A, xk, rk);
    subtract_v_v(b, rk, rk);

    // sk = rk
    copy_vector(rk, sk);

    double ak = 0; //  lenght of step
    double bk = 0; //  pom constant
    double rho = multiply_v_v(rk, rk);

    for (int i = 0; i < k; i++)
    {
        //  ak = rk'*rk / sk'*A*sk
        muptiply_m_v(A, sk, pom1);
        ak = rho / multiply_v_v(sk, pom1);
        
        //  rk_1 = rk - ak*A*sk
        // muptiply_m_v(A, sk, pom1);
        multiply_v_c(pom1, ak, pom1);
        subtract_v_v(rk, pom1, rk);
        
        //  xk_1 = xk + ak*sk
        multiply_v_c(sk, ak, pom1);
        addition_v_v(xk, pom1, xk);
        
        //  time to end?
        if (gen_vector_lenght(rk, A) < EPSILON)
            break;

        //  bk = rk_1' * rk_1 / rk' * rk
        bk = rho;
        rho = multiply_v_v(rk, rk);
        bk = rho/bk;
        // bk = rho / multiply_v_v(rk, rk);

        //  sk_1 = rk_1 + bk*sk
        multiply_v_c(sk, bk, sk);
        addition_v_v(rk, sk, sk);
        // rk = rk_1;
    }
}

//  main
int main(int argc, char *argv[])
{
    struct arguments arguments;

    //  default arguments
    arguments.format = FORMAT_2D; //  default 2d array atric format
    arguments.method = METHOD_GD; //  default gradient descent method
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
    printf("Matrix format: %s.\n", format_to_string(arguments.format));
    printf("Method: %s.\n", method_to_string(arguments.method));
    printf("------------------------------------------\n");
    printf("INPUTS READING:\n");
    printf("------------------------------------------\n");

    //  declare matrix
    double **matrix;

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
        if (read_matrix(matrix, matrix_input) == 0)
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

    //  read vector
    if (read_vector(vector, vector_input) == 0)
    {
        printf("filename %s does not exist", vector_input);
        return 1;
    }

    //  log inform
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
        gradient_descent_method(matrix, vector, output);
        break;
    case METHOD_CG:
        conjugate_gradient_method(matrix, vector, output);
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
    free(output);

    return 0;
}
