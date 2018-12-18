#include<stdio.h>
#include<math.h>
#define MAXN 5 /* 最大采样点个数 */
#define MAXM 10 /* 最大估算点个数 */

double get_li(double x, int i, double t[], int N)
{
    double xi = t[i];
    double result = 1;
    int j;
    for(j = 0; j < N; j++)
    {
        if(j == i)
            continue;
        double xj = t[j];
        result *= (x-xj)/(xi-xj);
    }
    return result;
}

double get_star1(double x, int i, double t[], int N)
{
    int j;
    double result = 0;
    for(j = 0; j<N; j++)
    {
        double xj = t[j];
        if(j == i || xj == x)
            continue;
        result += 1/(x-xj);
    }
    return result;
}

double get_star2(int i, double t[], int N)
{
    double result = 0;
    double xi = t[i];
    int j;
    for(j = 0; j < N; j++)
    {
        double xj = t[j];
        if(j == i || xi == xj)
            continue;
        result += 1/(xi-xj);
    }
    return result;
}


double get_hi_x(double x, int i, double t[], int N)
{
    double xi = t[i];
    double li_differential = get_star2(i, t, N);
    double li_square = pow(get_li(x, i, t, N),2);
    double hi_x = (1-2*li_differential*(x-xi)) * li_square;
    return hi_x;
}

double get_hi_x_hat(double x, int i, double t[], int N)
{
    double li_square = pow(get_li(x, i, t, N),2);
    double xi = t[i];
    double hi_x_hat = (x-xi)*li_square;
    return hi_x_hat;
}

double get_Hx(double x, double t[], double s[], double v[], int N)
{
    int i;
    double first = 0.0;
    double second = 0.0;
    for(i = 0; i < N; i++)
    {
        double yi = s[i];
        double yi_differential = v[i];
        first += yi * get_hi_x(x, i, t, N);
        second += yi_differential * get_hi_x_hat(x, i, t, N);
    }
    return first + second;
}


double get_Hx_differential(double x, double t[], double s[], double v[], int N)
{
    int i;
    double first = 0;
    double second = 0;
    double result = 0;
    for(i =0; i < N; i++)
    {
        double xi = t[i];
        double yi = s[i];
        double yi_differential = v[i];
        double li_square = pow(get_li(x, i, t, N),2);
        first = 2 * yi * li_square * (get_star1(x, i, t, N) * (1-2*get_star2(i, t, N)*(x-xi)) - get_star2(i, t, N));
        second = yi_differential * li_square * (1+2*get_star1(x, i, t, N)*(x-xi));
        result += first + second;
    }
    return result;
}
void Hermite_Interpolation( int N, double t[], double s[], double v[], int m, double ht[], double hs[], double hv[] )
{
    int i;
    for(i = 0; i < m; i++)
    {
        double x = ht[i];
        hs[i] = get_Hx(x, t, s, v, N);
        hv[i] = get_Hx_differential(x, t, s, v, N);
    }
}

int main()
{
  int N, m;
  double t[MAXN], s[MAXN], v[MAXN]; /* 用于构造 的数据 */
  double ht[MAXM], hs[MAXM], hv[MAXM]; /* 用 估算的数据 */
  int i;
  while ( scanf("%d", &N) != EOF ) {
    for ( i=0; i<N; i++ )
        scanf("%lf", &t[i]);
    for ( i=0; i<N; i++ )
        scanf("%lf", &s[i]);
    for ( i=0; i<N; i++ )
        scanf("%lf", &v[i]);
    scanf("%d", &m);
    for ( i=0; i<m; i++ )
      scanf("%lf", &ht[i]);
    Hermite_Interpolation( N, t, s, v, m, ht, hs, hv );
    for ( i=0; i<m; i++ )
      printf("%.4lf ", hs[i]);
    printf("\n");
    for ( i=0; i<m; i++ )
      printf("%4lf ", hv[i]);
    printf("\n\n");
  }
  return 0;
}

