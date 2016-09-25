#include <iostream>
#include <gmpxx.h>
#include <math.h>
#include <fstream>
#include <string>
using std::cout;
using std::endl;
int blum_mikali(mpz_class &T);
char* blum_mikali_byte(mpz_class &T, int N);
int bbs(mpz_class &r);
char* bbs_byte(mpz_class &r, int N);
double quantile(double a);
double chi_squared(double Z, int N);
bool criterion_equal_probability(char Y[], int m, double a);
bool criterion_sign_independency(char Y[], int m, double a);
bool criterion_homogeneity(char Y[], int m, double a, int r);
char* blum_mikali_array(mpz_class &T, int N); 
char* bbs_array(mpz_class &r, int N);
void check(char Y[], int m, double a, int r);
void left_shift(int L[], int size);
void right_shift(int L[], int size);
int built_in_rand();
int lehmer_low(int &L);
int lehmer_high(int &L);
int l20(int &x);
int l89(int x[3], int size);
int geffe(int &L1, int &L2, int &L3);
int wolfram(int &L);
int main() {
  mpz_class T,r;
  T = pow(2,12);
  r = pow(5,56);
  int N = 200000;
  double a = 0.01;
  int r_ch = 66;
  double alpa = 0.001;
  
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<unsigned int> uni(1, UINT32_MAX);
  
  char mas[131072] = "";
  int L;
  int L_mas[3];
  int L1, L2, L3;
 
  cout << "1) Built-in RNG" << endl;
  int seed = (time(nullptr)); // start seed
  srand(seed);
  for (auto i = 0; i < 131072; i++)
  {
    mas[i] = built_in_rand();
  }	
  check(mas, 131072, alpa, 24);
  
  cout << "2) Lemmer-Low" << endl;
  L = uni(rng) & -1;
  for (auto i = 0; i < 131072; i++)
  {
    mas[i] = lehmer_low(L);
  }
  check(mas, 131072, alpa, 24);
  
  cout << "3) Lemmer-High" << endl;
  L = uni(rng) & -1;
  for (auto i = 0; i < 131072; i++)
  {
    mas[i] = lehmer_high(L);
  }
  check(mas, 131072, alpa, 24);
  
  cout << "4) L20" << endl;
  L = uni(rng) & -1;
  for (auto i = 0; i < 131072; i++)
  {
    for (auto j = 0; j < 8; j++)
    {
      mas[i] <<= 1;
      mas[i] ^= l20(L);
    }
  }
  check(mas, 131072, alpa, 24);

  cout << "5) L89" << endl;
  for (auto i = 0; i < 3; i++)
  {
    L_mas[i] = 56312;
  }
  for (auto i = 0; i < 131072; i++)
  {
    for (auto j = 0; j < 8; j++)
    {
      mas[i] <<= 1;
      mas[i] ^= l89(L_mas, 3);
    }
  }
  check(mas, 131072, alpa, 24);
  
  cout << "6) Geffe" << endl;
  L1= uni(rng) & -1;
  L2= uni(rng) & -1;
  L3= uni(rng) & -1;
  for (auto i = 0; i < 131072; i++)
  {
    for (auto j = 0; j < 8; j++)
    {
      mas[i] <<= 1;
      mas[i] ^= geffe(L1, L2, L3);
    }
  }
  check(mas, 131072, alpa, 24);
  
  cout << "7) Wolfram" << endl;
  L = uni(rng) & -1;
  for (auto i = 0; i < 131072; i++)
  {
    for (auto j = 0; j < 8; j++)
    {
      mas[i] <<= 1;
      mas[i] ^= wolfram(L);
    }
  }
  check(mas, 131072, alpa, 24);  

  char* bmdata = blum_mikali_array(T,N);
  std::cout<<"passing test by blum_mikali"<<" coeficient a = "<<a<<std::endl;
  std::cout<<criterion_equal_probability(bmdata,N,a)<<std::endl;
  std::cout<<criterion_sign_independency(bmdata,N,a)<<std::endl;
  std::cout<<criterion_homogeneity(bmdata,N,a,r_ch)<<std::endl;
  std::cout<<std::endl;
  char* bbsdata = bbs_array(T,N);
  std::cout<<"passing test by bbs"<<" coeficient a = "<< a<<std::endl;
  std::cout<<criterion_equal_probability(bbsdata,N,a)<<std::endl;
  std::cout<<criterion_sign_independency(bbsdata,N,a)<<std::endl;
  std::cout<<criterion_homogeneity(bbsdata,N,a,r_ch)<<std::endl;
  std::cout<<std::endl;
  char* data_for_bm_byte = blum_mikali_byte(T,N);
  std::cout<<"passing test by blum_mikali_byte"<<" coeficient a = "<< a<<std::endl;
  std::cout<<criterion_equal_probability(data_for_bm_byte,N,a)<<std::endl;
  std::cout<<criterion_sign_independency(data_for_bm_byte,N,a)<<std::endl;
  std::cout<<criterion_homogeneity(data_for_bm_byte,N,a,r_ch)<<std::endl;
  std::cout<<std::endl;
  char* data_for_bbs_byte = bbs_byte(r,N);
  std::cout<<"passing test by bbs_byte"<<" coeficient a = "<< a<<std::endl;
  std::cout<<criterion_equal_probability(data_for_bbs_byte,N,a)<<std::endl; 
  std::cout<<criterion_sign_independency(data_for_bbs_byte,N,a)<<std::endl;
  std::cout<<criterion_homogeneity(data_for_bbs_byte,N,a,r_ch)<<std::endl;
  delete [] data_for_bbs_byte, data_for_bm_byte, bmdata, bbsdata;
  return 0;
}
void left_shift(int L[], int size)
{ 
  auto i = -1;
  while (++i < size-1)
  {
    L[i] <<= 1;
    L[i] ^= (L[i + 1] & 1) << 31;
  }
  L[size - 1] <<= 1;
}

void right_shift(int L[], int size)
{
  while (--size > 0)
  {
    L[size] >>= 1;
    L[size] ^= (L[size - 1] & 1) << 31;
  }
  L[0] >>= 1;
}

void check(char Y[], int m, double a, int r)
{
  criterion_equal_probability(Y, m, a);
  criterion_sign_independency(Y, m, a);
  criterion_homogeneity(Y, m, a, r);
}

int built_in_rand()
{
  return rand() % 256 + 1;
}

int lehmer_low(int &L)
{
  L = ((1 << 16) + 1) * L + 119;
  return  L & 255;
}

int lehmer_high(int &L)
{
  L = ((1 << 16) + 1) * L + 119;
  return  L >> 24;
}

int l20(int &L)
{
  auto out = L & 1;
  if (((L >> 16) ^ (L >> 14) ^ (L >> 10) ^ L) & 1)
  {
    L >>= 1;
    L ^= 1 << 19;
  }
  else
  {
    L >>= 1;
  }
  return out;
}

int l89(int L[3], int size)
{
  auto out = L[2] & 1;
  if (((L[1] >> 18) ^ L[2]) & 1)
  {
    right_shift(L, size);
    L[0] ^= 1 << 24;
  }
  else
  {
    right_shift(L, size);
  }
  return out;
}

int geffe(int &L1, int &L2, int &L3)
{
  auto x = L1 & 1;
  auto y = L2 & 1;
  auto s = L3 & 1;
  L1 >>= 1;
  L2 >>= 1;
  L3 >>= 1;
  L1 ^= ((L1 ^ L1 >> 2) & 1) << 11;
  L2 ^= ((L2 ^ L2 >> 1 ^ L2 >> 3 ^ L2 >> 4) & 1) << 9;
  L3 ^= ((L3 ^ L3 >> 3) & 1) << 10;
  return s*x ^ (1 ^ s)*y;
}

int wolfram(int &L)
{
  auto out = L & 1;
  L = (L << 1 ^ L >> 31) ^ (L | L >> 1 ^ (L & 1) << 31);
  return out;
}

char* blum_mikali_array(mpz_class &T, int N)
{
  char* out = new char[N];
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < 8; j++)
    {
      out[i] <<= 1;
      out[i] ^=blum_mikali(T);   
    }
  }
  return out;
}
int blum_mikali(mpz_class &T)
{
  //char* out_put = new char[N];
  mpz_class p("0xCEA42B987C44FA642D80AD9F51F10457690DEF10C83D0BC1BCEE12FC3B6093E3");
  mpz_class a("0x5B88C41246790891C095E2878880342E88C79974303BD0400B090FE38A688356");
  mpz_class one("1");
  mpz_class two("2");
  mpz_class c;
  c = (p-one)/2;
  mpz_powm(T.get_mpz_t(),a.get_mpz_t(),T.get_mpz_t(),p.get_mpz_t());  
 // for(int i = 0; i < N; i++)
 // {
    if(T < c)
    {
      //out_put[i] = 1;
     // mpz_powm(T.get_mpz_t(),a.get_mpz_t(),T.get_mpz_t(),p.get_mpz_t()); 
      return 1; 
    }
    else if(T >= c)
    { 
     // out_put[i] = 0;
     // mpz_powm(T.get_mpz_t(),a.get_mpz_t(),T.get_mpz_t(),p.get_mpz_t());
      return 0;
    }
  //}
  return -1;
}

char* blum_mikali_byte(mpz_class &T, int N)
{
  char* out_put = new char[N];
  mpz_class p("0xCEA42B987C44FA642D80AD9F51F10457690DEF10C83D0BC1BCEE12FC3B6093E3");
  mpz_class a("0x5B88C41246790891C095E2878880342E88C79974303BD0400B090FE38A688356");
  mpz_class one("1");
  mpz_class two("2");
  mpz_class c;
  c = (p-one)/2;
  mpz_class p_One;
  p_One = p-one;
  mpz_class b("256"), k;
 
  for(int i = 0; i < N; i++)
  {
    mpz_powm(T.get_mpz_t(),a.get_mpz_t(),T.get_mpz_t(),p.get_mpz_t());
    k  =  ((T * b)/(p-one));
    out_put[i] = k.get_ui();
  }
  return out_put;
}

char* bbs_array(mpz_class &r, int N)
{
  char* out = new char[N];
  int a =0;
  int b =0;
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < 8; j++)
    {
      out[i] <<= 1;
      out[i] ^= bbs(r);
    }
  }
  return out;
}

int bbs(mpz_class &r)
{
 // char* out_put = new char[N];
  mpz_class p("0xD5BBB96D30086EC484EBA3D7F9CAEB07");
  mpz_class q("0x425D2B9BFDB25B9CF6C416CC6E37B59C1F");
  mpz_class n = p*q;
  mpz_class x;
  mpz_class two("2");
 // for(int i = 0; i < N; i++)
 // {
    mpz_powm(r.get_mpz_t(),r.get_mpz_t(),two.get_mpz_t(),n.get_mpz_t());
    mpz_mod(x.get_mpz_t(),r.get_mpz_t(),two.get_mpz_t());
   // x = r.get_mpz_t()%two.get_mpz_t();
   // out_put[i] = x.get_ui();
 // }
  return x.get_ui();
}

char* bbs_byte(mpz_class &r,int N)
{
  char* out_put = new char[N];
  mpz_class p("0xD5BBB96D30086EC484EBA3D7F9CAEB07");
  mpz_class q("0x425D2B9BFDB25B9CF6C416CC6E37B59C1F");
  mpz_class n = p*q;
  mpz_class x;
  mpz_class two("2"); 
  mpz_class b("256");
 
  for(int i = 0; i < N; i++)
  {
    mpz_powm(r.get_mpz_t(),r.get_mpz_t(),two.get_mpz_t(),n.get_mpz_t());
    //x = r.get_ui()%b.get_ui();
    mpz_mod(x.get_mpz_t(),r.get_mpz_t(),b.get_mpz_t());
    out_put[i] = x.get_ui();
  }
  return out_put;
}
/************************************************************************
*********************************TESTS***********************************
*************************************************************************/
double quantile(double a)
{
  double mas_a[8] = { 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999};
  double mas_q[8] = { 0.841621234, 1.036433389, 1.281551566, 1.644853627, 2.326347874, 3.090232306, 3.719016485, 4.264890794 };
  int count = 0;
  while (mas_a[count] < 1 - a) count++;
  return mas_q[count];
}

double chi_squared(double Z, int N)
{
  return sqrt(2 * (N - 1))*Z + N - 1;
} 
bool criterion_equal_probability(char Y[], int m, double a)
{
  int v[256] = {};
  double chi = 0;
  auto n = m / 256;
  for (auto i = 0; i < m; i++)
  {
    v[Y[i] & 255]++;
  }
  for (auto i = 0; i < 256; i++)
  {
    chi += pow(v[i] - n, 2) / n;
  }
  chi = int(chi);
  int standart = chi_squared(quantile(a), 256);
  std::cout << "a)Equal probability criterion:" << std::endl;
  std::cout << "Real: " << chi << "\n Standart: " << standart <<std:: endl;
  std::cout << "Verdict: ";
  if (chi <= standart && chi>= 0 )
    std::cout << "true" <<std::endl;
  else
    std::cout << "false" <<std::endl;
  return (chi <= standart && chi >= 0);
}

bool criterion_sign_independency(char Y[], int m, double a)
{
  int v[256][512] = {};
  int mas_v[256] = {};
  int mas_a[256] = {};
  double chi = 0;
  auto n = m / 2;
  for (auto i = 0; i < n; i++)
  {
    v[Y[2 * i - 1] & 255][Y[2 * i] & 255]++;
    mas_v[Y[2 * i - 1] & 255]++;
    mas_a[Y[2 * i] & 255]++;
  }
  for (auto i = 0; i < 256; i++)
  {
    for (auto j = 0; j < 256; j++)
    {
      chi += pow(v[i][j], 2) / (mas_v[i] * mas_a[j]);
    }
  }
  chi--;
  chi *= n;
  chi = int(chi);
  int standart = chi_squared(quantile(a), 255 * 255 + 1);
  std::cout << "b)Sign independency criterion:" <<std::endl;
  std::cout << "Real: " << chi << "\n Standart: " << standart <<std::endl;
  std::cout << "Verdict: ";
  if (chi <= standart && chi >= 0 )
    std::cout << "true" <<std::endl;
  else
    std::cout << "false" << std::endl;
  return (chi <= standart && chi >= 0);
}

bool criterion_homogeneity(char Y[], int m, double a, int r)
{
  double chi = 0;
  auto m_ = m / r;
  auto n = m_*r;
  int v[256][512] = {};
  int mas_v[256] = {};
  for (auto i = 0; i < n; i++)
  {
    v[Y[i] & 255][i / m_]++;
    mas_v[Y[i] & 255]++;
  }
  for (auto i = 0; i < 256; i++)
  {
    for (auto j = 0; j < r; j++)
    {
      chi += pow(v[i][j], 2) / (mas_v[i] * m_);
    }
  }
  chi--;
  chi *= n;
  chi = int(chi);
  int standart = chi_squared(quantile(a), 255 * (r - 1) + 1);
  std::cout << "c)Binary homogeneity criterion:" <<std::endl;
  std::cout << "Real: " << chi << "\nStandart: " << standart <<std::endl;
  std::cout << "Verdict: ";
  if (chi <= standart && chi >= 0)
    std::cout << "true" <<std::endl;
  else
    std::cout << "false" << std::endl;
  return (chi <= standart && chi >= 0) ;
}



