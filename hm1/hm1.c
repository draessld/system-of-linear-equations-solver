#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <argp.h>
#include <stdbool.h>

#define EPSILON 1e-5

size_t n;

//  parsing setting
const char *argp_program_version = "1.0";
const char *argp_program_bug_address = "<draesdom@cvut.cz>";
static char doc[] = "Solve algebraic equations";
static char args_doc[] = "[FILENAMES]...";
static struct argp_option options[] = {
    {"2da", 'a', 0, 0, "Save matrix as 2d array format."},
    {"gdm", 'g', 0, 0, "Use gradient descent method."},
    {"output", 'o', "OUTFILE", 0, "Output file path."},
    {"size", 'n', "N", 0, "squared matrix size"},
    {0}};

struct arguments
{
    enum method
    {
        METHOD_GD
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
    //  output files
    case 'o':
        arguments->outfile = arg;
        break;
    //  other
    case 'n':
        arguments->count_n = false;
        n = atoi(arg);
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
void multiply_v_c(double *a, double c)
{
    for (int i = 0; i < n; i++)
        a[i] *= c;
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

//  gradient descent method
void gradient_descent_method(double **matrix, double *b, double *x)
{

    //  initialization vector x
    double pom[n];
    for (int i = 0; i < n; i++)
    {
        x[i] = 0;
        pom[i] = 0;
    }

    double *gk = (double *)malloc(n * sizeof(double));
    muptiply_m_v(matrix, x, pom);
    subtract_v_v(b, pom, pom);

    copy_vector(pom, gk);
    double lgo = vector_lenght(gk);
    double a = 0; //  lenght of step
    double e = 0; //  difference

    //  lets do it! :))
    do
    {
        double y = multiply_v_v(gk, gk);
        // printf("y: %lf \n",y);
        muptiply_m_v(matrix, gk, pom);
        // print_v(pom);
        a = y / multiply_v_v(gk, pom);
        // printf("a: %lf \n", a);
        multiply_v_c(gk, a);
        // print_v(gk);
        addition_v_v(x, gk, x);
        muptiply_m_v(matrix, x, pom);
        subtract_v_v(b, pom, gk);
        // print_v(x);
        double lgk = vector_lenght(gk);
        if (lgo == 0)
            e = lgk;
        else
            e = lgk / lgo;
        // printf("e: %lf \n",e);
    } while (e > EPSILON);

    free(gk);
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

    //  switch by method
    switch (arguments.method)
    {
    case METHOD_GD:
        gradient_descent_method(matrix, vector, output);
        break;

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
