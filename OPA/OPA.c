#include<stdio.h>
#include<math.h>

#define MAXM 200 /* 采样点个数上限*/
#define MAXN 5 /* 拟合多项式阶数上限*/

/* 这里仅给出2个测试用例，正式裁判程序可能包含更多测试用例 */

double f1(double x)
{
  return sin(x);
}

double f2(double x)
{
  return exp(x);
}

int OPA( double (*f)(double t), int m, double x[], double w[], double c[], double *eps )
{
    //--------------definition----------------
    double err;
    double phi0[MAXN+1],phi1[MAXN+1],phi2[MAXN+1];
    double alpha, beta;
    double a[MAXN+1];
    for(int i=0; i<=MAXN;i++)
    {
        phi0[i] = 0;
        phi1[i] = 0;
        phi2[i] = 0;
    }

    //--------------Step 1------------------------//

    phi0[0] = 1;
    double y_y=0, phi0_y=0, phi0_phi0=0;
    for(int i=0;i<m;i++)
    {
        phi0_y += w[i] * f(x[i]);
        phi0_phi0 += w[i];
    }
    a[0] = phi0_y / phi0_phi0;
    for(int i=0;i<m;i++)
        y_y += w[i] * f(x[i]) * f(x[i]);
    err = y_y - a[0] * phi0_y;
    for(int i=0; i<=MAXN; i++)
        c[i] = a[0] * phi0[i];
    printf("err:%f\n",err);
    //---------------Step 2------------------------//

    double xphi0_phi0=0, phi1_y=0, phi1_phi1=0;
    phi0_phi0 = 0;
    for(int i=0; i<m; i++)
    {
        xphi0_phi0 += w[i] * x[i];
        phi0_phi0 += w[i];
    }
    alpha = xphi0_phi0 / phi0_phi0;
    phi1[0] = -alpha;
    phi1[1] = 1;
    for(int i=0; i<m; i++)
    {
        phi1_y += w[i] * f(x[i]) * (x[i]-alpha);
        phi1_phi1 += w[i]* pow((x[i]-alpha),2);
    }
    a[1] = phi1_y / phi1_phi1;
    err -= a[1] * phi1_y;
    for(int i=0; i<=MAXN; i++)
        c[i] += a[1] * phi1[i];
    printf("err:%f\n",err);
    //--------------Step 345678-------------------//

    int k = 1;
    while( (k<MAXN) && (fabs(err)>=(*eps)))
    {
        printf("err:%f\n",err);
        k++;
        double xphi1_phi1=0;
        phi1_phi1 = 0, phi0_phi0 = 0;
        for(int i=0; i<m; i++)
        {
            double Phi0 = 0, Phi1 = 0;
            for(int j=0; j<=MAXN; j++)
            {
                Phi0 += phi0[j] * pow(x[i],j);
                Phi1 += phi1[j] * pow(x[i],j);
            }
            phi1_phi1 += w[i] * Phi1 * Phi1;
            xphi1_phi1 += w[i] * x[i] * Phi1 * Phi1;
            phi0_phi0 += w[i] * Phi0 * Phi0;
        }
        alpha = xphi1_phi1 / phi1_phi1;
        beta = phi1_phi1 / phi0_phi0;
        printf("beta:%f\n",beta);
        phi2[0] = -alpha*phi1[0] - beta*phi0[0];
        for(int i=1; i<=MAXN; i++)
            phi2[i] = phi1[i-1] - alpha*phi1[i] - beta*phi0[i];

        double phi2_y=0, phi2_phi2=0;
        for(int i=0;i<m;i++)
        {
            double Phi2=0;
            for(int j=0; j<=MAXN; j++)
                Phi2 += phi2[j] * pow(x[i],j);
            phi2_y += w[i]* Phi2 * f(x[i]);
            phi2_phi2 += w[i] * Phi2 * Phi2;
        }
        a[k] = phi2_y / phi2_phi2;
        err -= a[k] * phi2_y;
        for(int i=0; i<=MAXN; i++)
        {
            c[i] += a[k] * phi2[i];
            phi0[i] = phi1[i];
            phi1[i] = phi2[i];
        }
    }
    *eps = err;
    return k;
}
void print_results( int n, double c[], double eps)
{ /* 用于输出结果 */
  int i;

  printf("%d\n", n);
  for (i=0; i<=n; i++)
    printf("%8.4e ", c[i]);
  printf("\n");
  printf("error = %6.2e\n", eps);
  printf("\n");
}

int main()
{
  int m, i, n;
  double x[MAXM], w[MAXM], c[MAXN+1], eps;

  /* 测试1 */
  m = 90;
  for (i=0; i<m; i++) {
    x[i] = 3.1415926535897932 * (double)(i+1) / 180.0;
    w[i] = 1.0;
  }
  eps = 0.001;
  n = OPA(f1, m, x, w, c, &eps);
  print_results(n, c, eps);

  /* 测试2 */
  m = 200;
  for (i=0; i<m; i++) {
    x[i] = 0.01*(double)i;
    w[i] = 1.0;
  }
  eps = 0.001;
  n = OPA(f2, m, x, w, c, &eps);
  print_results(n, c, eps);

  return 0;
}

