#include<stdio.h>
#include<math.h>
#include<time.h>
double f0( double x, double l, double t )
{ /* 弧长积分的核函数 */
  return sqrt(1.0+l*l*t*t*cos(t*x)*cos(t*x));
}

double T(int n, double a, double b, double (*f)(double x, double y, double z), double l, double t)
{
    double T0 = 0.5*(f(a,l,t)+f(b,l,t));
    if(n==0)
        return T0;
    double h = (b-a)/2;
    double T1 = 0.5*T0+0.5*f(a+h,l,t);
    if(n==1)
        return T1;
    else
    {
        double Tm;
        int m = 1;
        while(m<n)
        {
            m = m+1;
            h = h/2;
            double sum = 0;
            sum += f(a+h,l,t);
            for(int i=2;i<=pow(2,m-1);i++)
                sum += f(a+h+(i-1)*2*h,l,t);
            Tm = 0.5*T1+ h*sum;
            T1 = Tm;
        }
        return Tm;
    }
}

double get_temp(int m, int k, double a, double b, double (*f)(double x, double y, double z), double l, double t)
{
    if(m==0)
        return T(k, a, b, f, l, t);
    else return pow(4,m)*get_temp(m-1,k+1,a,b,f,l,t)/(pow(4,m)-1) - get_temp(m-1,k,a,b,f,l,t)/(pow(4,m)-1);
}

double Integral(double a, double b, double (*f)(double x, double y, double z), double TOL, double l, double t)
{
    int m = 0;
    int k = 0;
    double temp1 = get_temp(m,k,a,b,f,l,t);
    double temp2 = get_temp(m+1,k,a,b,f,l,t);
    while(fabs(temp2-temp1)>=TOL)
    {
        m += 1;
        temp1 = temp2;
        temp2 = get_temp(m+1,k,a,b,f,l,t);
    }
    return 0.01*temp2;
}

int main()
{
  double a=0.0, b, TOL=0.005, l, t;
  while (scanf("%lf %lf %lf", &l, &b, &t) != EOF)
  {
    clock_t start, end;
    start = clock();
    double result = Integral(a, b, f0, TOL, l, t);
    end = clock();
    printf("%.2f\n", result);
    printf("time:%f\n",(double)(end-start)/CLOCKS_PER_SEC);
  }
  return 0;
}
