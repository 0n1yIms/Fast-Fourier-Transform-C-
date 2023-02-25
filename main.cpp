#include <Halide.h>
#include <iostream>

#define PI 3.1415926535897932384626433832795f

using namespace std;
using namespace Halide;

struct Cmpx
{
  float real;
  float imag;
  Cmpx(float r)
  {
    real = r;
    imag = 0.f;
  }
  Cmpx(float r, float i)
  {
    real = r;
    imag = i;
  }
  Cmpx()
  {
    real = 0;
    imag = 0;
  }

};


#define N 16
#define M 4
#define LVLS 5
/**
 * N = 16
 * m = 4
 * lvls = 5
 * nodes = [16, 8, 4, 2, 1]
 * fqLens = [1, 2, 4, 8, 16]
 * */

// int lvl0Idx[] = {0, 4, 2, 6, 1, 5, 3, 7};      N=8
// int nodes[] = {8, 4, 2, 1};
// int fqLens[] = {1, 2, 4, 8};

int lvl0Idx[] = {0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};
int nodes[] = {16, 8, 4, 2, 1};
int fqLens[] = {1, 2, 4, 8, 16};

int idxFt(int lvl, int fq, int node){
  return lvl * N + fq * nodes[lvl] + node;
}

Cmpx* fftAlg(float *img) {
  Cmpx *lvlF = new Cmpx[N * LVLS];
  
  for (int i = 0; i < N; i++)
    lvlF[idxFt(0, 0, i)] = img[lvl0Idx[i]];

  for (int lvl = 1; lvl < LVLS; lvl++) {
    int node = nodes[lvl];
    int fqLen = fqLens[lvl];

    for (int fq = 0; fq < fqLen; fq++)
    {
      for (int nd = 0; nd < node; nd++)
      {
        Cmpx a = lvlF[idxFt(lvl - 1,fq % fqLens[lvl - 1], 2 * nd)];
        Cmpx b = lvlF[idxFt(lvl - 1,fq % fqLens[lvl - 1], 2 * nd + 1)];

        Cmpx e = {cos(2* PI * fq / fqLen), -1 * sin(2 * PI * fq / fqLen)};
        Cmpx mult = {b.real * e.real - b.imag * e.imag, b.real * e.imag + b.imag * e.real};

        lvlF[idxFt(lvl,fq, nd)] = {a.real + mult.real, a.imag + mult.imag};
      }      
    }
  }
  
  lvlF = &lvlF[idxFt(M, 0, 0)];

  return lvlF;
}

float* ifftAlg(Cmpx *ft)
{
  for (int i = 0; i < N; i++)
    ft[i] = {ft[i].real, -1 * ft[i].imag};

  Cmpx lvlF[LVLS * N];
  
  for (int i = 0; i < N; i++)
    lvlF[idxFt(0, 0, i)] = ft[lvl0Idx[i]];

  for (int lvl = 1; lvl < LVLS; lvl++) {
    int node = nodes[lvl];
    int fqLen = fqLens[lvl];

    for (int fq = 0; fq < fqLen; fq++)
    {
      for (int nd = 0; nd < node; nd++)
      {
        Cmpx a = lvlF[idxFt(lvl - 1,fq % fqLens[lvl - 1], 2 * nd)];
        Cmpx b = lvlF[idxFt(lvl - 1,fq % fqLens[lvl - 1], 2 * nd + 1)];

        Cmpx e = {cos(2* PI * fq / fqLen), -1 * sin(2 * PI * fq / fqLen)};
        Cmpx mult = {b.real * e.real - b.imag * e.imag, b.real * e.imag + b.imag * e.real};

        lvlF[idxFt(lvl,fq, nd)] = {a.real + mult.real, a.imag + mult.imag};
      }      
    }
  }
  float *ifft = new float[N];
  for (int  i = 0; i < N; i++)
  {
    ifft[i] = lvlF[idxFt(M, i, 0)].real / float(N);
  }
  
  return ifft;
}



/*
0: 0 0
1: 8 0
2: 4 0
3: 12 0
4: 2 0
5: 10 0
6: 6 0
7: 14 0
8: 1 0
9: 9 0
10: 5 0
11: 13 0
12: 3 0
13: 11 0
14: 7 0
15: 15 0

16: 8 0
17: 16 0
18: 12 0
19: 20 0
20: 10 0
21: 18 0
22: 14 0
23: 22 0
24: -8 6.99382e-07
25: -8 1.04907e-06
26: -8 8.74228e-07
27: -8 1.22392e-06
28: -8 7.86805e-07
29: -8 1.1365e-06
30: -8 9.61651e-07
31: -8 1.31134e-06

32: 24 0
33: 32 0
34: 28 0
35: 36 0
36: -8 8
37: -8 8
38: -8 8
39: -8 8
40: -8 1.39876e-06
41: -8 1.74846e-06
42: -8 1.57361e-06
43: -8 1.9233e-06
44: -8 -8
45: -8 -8
46: -8 -8
47: -8 -8

48: 56 0
49: 64 0
50: -8 19.3137
51: -8 19.3137
52: -8 8
53: -8 8
54: -8 3.31371
55: -8 3.31371
56: -8 2.79753e-06
57: -8 3.14722e-06
58: -8 -3.31371
59: -8 -3.31371
60: -8 -8
61: -8 -8
62: -8.00001 -19.3137
63: -8.00001 -19.3137

64: 120 0
65: -7.99999 40.2187
66: -7.99999 19.3137
67: -8 11.9728
68: -8 8
69: -8 5.34543
70: -8 3.31371
71: -8 1.5913
72: -8 5.59506e-06
73: -8 -1.5913
74: -8 -3.31371
75: -8 -5.34543
76: -8 -8
77: -8.00001 -11.9728
78: -8.00001 -19.3137
79: -8.00001 -40.2187

*/

Func fftAlg(Func img) {
  Var x, y, z;
  Func nodes, fqLens;
  nodes(x) = select(x == 0, 16,
                   x == 1, 8,
                   x == 2, 4,
                   x == 3, 2,
                   x == 4, 1, 0);
  fqLens(x) = select(x == 0, 1,
                   x == 1, 2,
                   x == 2, 4,
                   x == 3, 8,
                   x == 4, 16, 0);


  Func lvl0("lvl0"), lvl1("lvl1"), lvl2("lvl2"), lvl3("lvl3"), lvl4("lvl4"), out("out");
  {
    Var x;
    //0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
    lvl0(x) = img(select(x == 0, 0,
                    x == 1, 8,
                    x == 2, 4,
                    x == 3, 12,
                    x == 4, 2,
                    x == 5, 10,
                    x == 6, 6,
                    x == 7, 14,
                    x == 8, 1,
                    x == 9, 9,
                    x == 10, 5,
                    x == 11, 13,
                    x == 12, 3,
                    x == 13, 11,
                    x == 14, 7,
                    x == 15, 15, 0));
    lvl0.compute_root();
    // lvl0.trace_stores();
  }
  //lvl1
  //nodes {16, 8, 4, 2, 1};
  //fqLens {1, 2, 4, 8, 16};
  {
    int node = 8;
    int fqLen = 2;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Expr a = lvl0(2 * nd);
    Expr b = lvl0(2 * nd + 1);

    // Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Expr e = cos(2.f * PI * fq / fqLen);
    Expr mult = b * e;

    lvl1(fq, nd) = a + mult;
    lvl1.compute_root();
    lvl1.vectorize(nd, 8);
    // lvl1.trace_stores();
    // lvl1.realize({2, 8});
  }
  //lvl2
  {
    int node = 4;
    int fqLen = 4;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Expr a = lvl1(fq % fqLens(1), 2 * nd);
    Expr b = lvl1(fq % fqLens(1), 2 * nd + 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b * e[0], b * e[1]};
    // Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl2(fq, nd) = Tuple(a + mult[0], mult[1]);
    lvl2.compute_root();
    // lvl2.trace_stores();
    // lvl2.realize({4, 4});
  }
  //lvl3
  {
    int node = 2;
    int fqLen = 8;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Tuple a = lvl2(fq % fqLens(2), 2 * nd);
    Tuple b = lvl2(fq % fqLens(2), 2 * nd + 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl3(fq, nd) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl3.compute_root();
    // lvl3.trace_stores();
  }
  //lvl4
  {
    int node = 1;
    int fqLen = 16;
    // RDom fq(0, fqLen);
    Var fq;

    Tuple a = lvl3(fq % fqLens(3), 0);
    Tuple b = lvl3(fq % fqLens(3), 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl4(fq) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl4.compute_root();
    cout << endl << "print lvl4" << endl;
    lvl4.trace_stores();
  }
  out(y, x) = select(y == 0, lvl4(x)[0], lvl4(x)[1]);

  return out;
}


Func ifftAlg(Func img) {
  Var x, y, z;

  Func nodes, fqLens;
  nodes(x) = select(x == 0, 16,
                   x == 1, 8,
                   x == 2, 4,
                   x == 3, 2,
                   x == 4, 1, 0);
  fqLens(x) = select(x == 0, 1,
                   x == 1, 2,
                   x == 2, 4,
                   x == 3, 8,
                   x == 4, 16, 0);


  Func lvl0("lvl0"), lvl1("lvl1"), lvl2("lvl2"), lvl3("lvl3"), lvl4("lvl4"), out("out");
  {
    Var x;
    //0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
    Expr sel = select(x == 0, 0,
                    x == 1, 8,
                    x == 2, 4,
                    x == 3, 12,
                    x == 4, 2,
                    x == 5, 10,
                    x == 6, 6,
                    x == 7, 14,
                    x == 8, 1,
                    x == 9, 9,
                    x == 10, 5,
                    x == 11, 13,
                    x == 12, 3,
                    x == 13, 11,
                    x == 14, 7,
                    x == 15, 15, 0);
    lvl0(x) = img(sel);
    lvl0.compute_root();
    // lvl0.trace_stores();
  }
  //lvl1
  //nodes {16, 8, 4, 2, 1};
  //fqLens {1, 2, 4, 8, 16};
  {
    int node = 8;
    int fqLen = 2;
    Var fq, nd;

    Tuple a = lvl0(2 * nd);
    Tuple b = lvl0(2 * nd + 1);

    // Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl1(fq, nd) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl1.compute_root();
    lvl1.vectorize(nd, 8);
    // lvl1.trace_stores();
    // lvl1.realize({2, 8});
  }
  //lvl2
  {
    int node = 4;
    int fqLen = 4;
    Var fq, nd;

    Tuple a = lvl1(fq % fqLens(1), 2 * nd);
    Tuple b = lvl1(fq % fqLens(1), 2 * nd + 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    // Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl2(fq, nd) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl2.compute_root();
    // lvl2.trace_stores();
    // lvl2.realize({4, 4});
  }
  //lvl3
  {
    int node = 2;
    int fqLen = 8;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Tuple a = lvl2(fq % fqLens(2), 2 * nd);
    Tuple b = lvl2(fq % fqLens(2), 2 * nd + 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    // Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl3(fq, nd) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl3.compute_root();
    // lvl3.trace_stores();
  }
  //lvl4
  {
    int node = 1;
    int fqLen = 16;
    // RDom fq(0, fqLen);
    Var fq;

    Tuple a = lvl3(fq % fqLens(3), 0);
    Tuple b = lvl3(fq % fqLens(3), 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    // Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl4(fq) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl4.compute_root();
    cout << endl << "print lvl4" << endl;
    lvl4.trace_stores();
  }
  out(x) = lvl4(x)[0] / 16.f;

  return out;
}


void ftN()
{
  float fun[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

  Cmpx *fft = fftAlg(&fun[0]);

  cout << "ft: " << endl;
  for (int i = 0; i < N; i++)
  {
    cout << fft[i].real << "   +    ";
    cout << fft[i].imag << "i" << endl;
  }

  float *ifft = ifftAlg(fft);

  cout << "ift: " << endl;
  for (int i = 0; i < N; i++)
    cout << ifft[i] << endl;
}
void hft()
{
  Buffer<float> funB(16);
  for (int i = 0; i < 16; i++)
  {
    funB(i) = float(i);
  }
  Func fun;
  Var x, y;
  fun(x) = funB(x);

  Buffer<float> fft = fftAlg(fun).realize({2, 16});
  

  cout << "ft: " << endl;
  for (int i = 0; i < N; i++)
  {
    cout << fft(0, i) << "   +    ";
    cout << fft(1, i) << "i" << endl;
    // cout << fft(i).real << "   +    ";
    // cout << fft(i).imag << "i" << endl;
  }

  Func ftF;
  ftF(x) = Tuple(fft(0, x), fft(1, x));
  Buffer<float> ifft = ifftAlg(ftF).realize({16});

  cout << "ift: " << endl;
  for (int i = 0; i < N; i++)
    cout << ifft(i) << endl;
}



Func fftAlg2D(Func img) {
  Var x, y, z;
  Func nodes, fqLens;
  nodes(x) = select(x == 0, 16,
                   x == 1, 8,
                   x == 2, 4,
                   x == 3, 2,
                   x == 4, 1, 0);
  fqLens(x) = select(x == 0, 1,
                   x == 1, 2,
                   x == 2, 4,
                   x == 3, 8,
                   x == 4, 16, 0);


  Func lvl0("lvl0"), lvl1("lvl1"), lvl2("lvl2"), lvl3("lvl3"), lvl4("lvl4"), out("out");
  {
    //0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
    Expr sel = select(x == 0, 0,
                    x == 1, 8,
                    x == 2, 4,
                    x == 3, 12,
                    x == 4, 2,
                    x == 5, 10,
                    x == 6, 6,
                    x == 7, 14,
                    x == 8, 1,
                    x == 9, 9,
                    x == 10, 5,
                    x == 11, 13,
                    x == 12, 3,
                    x == 13, 11,
                    x == 14, 7,
                    x == 15, 15, 0);
    lvl0(x, y) = img(sel, y);
    lvl0.compute_root();
    // lvl0.trace_stores();
  }
  //lvl1
  //nodes {16, 8, 4, 2, 1};
  //fqLens {1, 2, 4, 8, 16};
  {
    int node = 8;
    int fqLen = 2;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Expr a = lvl0(2 * nd, y);
    Expr b = lvl0(2 * nd + 1, y);

    // Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Expr e = cos(2.f * PI * fq / fqLen);
    Expr mult = b * e;

    lvl1(fq, nd, y) = a + mult;
    lvl1.compute_root();
    lvl1.vectorize(nd, 8);
    // lvl1.trace_stores();
    // lvl1.realize({2, 8});
  }
  //lvl2
  {
    int node = 4;
    int fqLen = 4;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Expr a = lvl1(fq % fqLens(1), 2 * nd, y);
    Expr b = lvl1(fq % fqLens(1), 2 * nd + 1, y);

    Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b * e[0], b * e[1]};
    // Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl2(fq, nd, y) = Tuple(a + mult[0], mult[1]);
    lvl2.compute_root();
    // lvl2.trace_stores();
    // lvl2.realize({4, 4});
  }
  //lvl3
  {
    int node = 2;
    int fqLen = 8;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Tuple a = lvl2(fq % fqLens(2), 2 * nd, y);
    Tuple b = lvl2(fq % fqLens(2), 2 * nd + 1, y);

    Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl3(fq, nd, y) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl3.compute_root();
    // lvl3.trace_stores();
  }
  //lvl4
  {
    int node = 1;
    int fqLen = 16;
    // RDom fq(0, fqLen);
    Var fq;

    Tuple a = lvl3(fq % fqLens(3), 0, y);
    Tuple b = lvl3(fq % fqLens(3), 1, y);

    Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl4(fq, y) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl4.compute_root();
    cout << endl << "print lvl4" << endl;
    lvl4.trace_stores();
  }

  Func lvl0_2("lvl0_2"), lvl1_2("lvl1_2"), lvl2_2("lvl2_2"), lvl3_2("lvl3_2"), lvl4_2("lvl4_2"), out_2("out_2");
  {
    //0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
    Expr sel = select(y == 0, 0,
                    y == 1, 8,
                    y == 2, 4,
                    y == 3, 12,
                    y == 4, 2,
                    y == 5, 10,
                    y == 6, 6,
                    y == 7, 14,
                    y == 8, 1,
                    y == 9, 9,
                    y == 10, 5,
                    y == 11, 13,
                    y == 12, 3,
                    y == 13, 11,
                    y == 14, 7,
                    y == 15, 15, 0);
    lvl0_2(x, y) = lvl4(x, sel);
    lvl0_2.compute_root();
    // lvl0.trace_stores();
  }
  //lvl1
  //nodes {16, 8, 4, 2, 1};
  //fqLens {1, 2, 4, 8, 16};
  {
    int node = 8;
    int fqLen = 2;
    Var fq, nd;

    Tuple a = lvl0_2(x, 2 * nd);
    Tuple b = lvl0_2(x, 2 * nd + 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl1_2(x, fq, nd) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl1_2.compute_root();
    lvl1_2.vectorize(nd, 8);
    // lvl1.trace_stores();
    // lvl1.realize({2, 8});
  }
  //lvl2
  {
    int node = 4;
    int fqLen = 4;
    Var fq, nd;

    Tuple a = lvl1_2(x, fq % fqLens(1), 2 * nd);
    Tuple b = lvl1_2(x, fq % fqLens(1), 2 * nd + 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    // Tuple mult = {b * e[0], b * e[1]};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl2_2(x, fq, nd) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl2_2.compute_root();
    // lvl2.trace_stores();
    // lvl2.realize({4, 4});
  }
  //lvl3
  {
    int node = 2;
    int fqLen = 8;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Tuple a = lvl2_2(x, fq % fqLens(2), 2 * nd);
    Tuple b = lvl2_2(x, fq % fqLens(2), 2 * nd + 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};


    lvl3_2(x, fq, nd) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl3_2.compute_root();
    // lvl3.trace_stores();
  }
  //lvl4
  {
    int node = 1;
    int fqLen = 16;
    // RDom fq(0, fqLen);
    Var fq;

    Tuple a = lvl3_2(x, fq % fqLens(3), 0);
    Tuple b = lvl3_2(x, fq % fqLens(3), 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), -1.f * sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl4_2(x, fq) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl4_2.compute_root();
    cout << endl << "print lvl4" << endl;
    lvl4_2.trace_stores();
  }
  out(x, y) = lvl4_2(x, y);

  return out;
}



Func ifftAlg2D(Func img) {
  Var x, y, z;
  Func nodes, fqLens;
  nodes(x) = select(x == 0, 16,
                   x == 1, 8,
                   x == 2, 4,
                   x == 3, 2,
                   x == 4, 1, 0);
  fqLens(x) = select(x == 0, 1,
                   x == 1, 2,
                   x == 2, 4,
                   x == 3, 8,
                   x == 4, 16, 0);


  Func lvl0("lvl0"), lvl1("lvl1"), lvl2("lvl2"), lvl3("lvl3"), lvl4("lvl4"), out("out");
  {
    //0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
    Expr sel = select(x == 0, 0,
                    x == 1, 8,
                    x == 2, 4,
                    x == 3, 12,
                    x == 4, 2,
                    x == 5, 10,
                    x == 6, 6,
                    x == 7, 14,
                    x == 8, 1,
                    x == 9, 9,
                    x == 10, 5,
                    x == 11, 13,
                    x == 12, 3,
                    x == 13, 11,
                    x == 14, 7,
                    x == 15, 15, 0);
    lvl0(x, y) = img(sel, y);
    lvl0.compute_root();
    // lvl0.trace_stores();
  }
  //lvl1
  //nodes {16, 8, 4, 2, 1};
  //fqLens {1, 2, 4, 8, 16};
  {
    int node = 8;
    int fqLen = 2;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Tuple a = lvl0(2 * nd, y);
    Tuple b = lvl0(2 * nd + 1, y);

    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};


    lvl1(fq, nd, y) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl1.compute_root();
    lvl1.vectorize(nd, 8);
    // lvl1.trace_stores();
    // lvl1.realize({2, 8});
  }
  //lvl2
  {
    int node = 4;
    int fqLen = 4;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Tuple a = lvl1(fq % fqLens(1), 2 * nd, y);
    Tuple b = lvl1(fq % fqLens(1), 2 * nd + 1, y);

    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};


    lvl2(fq, nd, y) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl2.compute_root();
    // lvl2.trace_stores();
    // lvl2.realize({4, 4});
  }
  //lvl3
  {
    int node = 2;
    int fqLen = 8;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Tuple a = lvl2(fq % fqLens(2), 2 * nd, y);
    Tuple b = lvl2(fq % fqLens(2), 2 * nd + 1, y);

    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl3(fq, nd, y) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl3.compute_root();
    // lvl3.trace_stores();
  }
  //lvl4
  {
    int node = 1;
    int fqLen = 16;
    // RDom fq(0, fqLen);
    Var fq;

    Tuple a = lvl3(fq % fqLens(3), 0, y);
    Tuple b = lvl3(fq % fqLens(3), 1, y);

    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl4(fq, y) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl4.compute_root();
    cout << endl << "print lvl4" << endl;
    lvl4.trace_stores();
  }

  Func lvl0_2("lvl0_2"), lvl1_2("lvl1_2"), lvl2_2("lvl2_2"), lvl3_2("lvl3_2"), lvl4_2("lvl4_2"), out_2("out_2");
  {
    //0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
    Expr sel = select(y == 0, 0,
                    y == 1, 8,
                    y == 2, 4,
                    y == 3, 12,
                    y == 4, 2,
                    y == 5, 10,
                    y == 6, 6,
                    y == 7, 14,
                    y == 8, 1,
                    y == 9, 9,
                    y == 10, 5,
                    y == 11, 13,
                    y == 12, 3,
                    y == 13, 11,
                    y == 14, 7,
                    y == 15, 15, 0);
    lvl0_2(x, y) = lvl4(x, sel);
    lvl0_2.compute_root();
    // lvl0.trace_stores();
  }
  //lvl1
  //nodes {16, 8, 4, 2, 1};
  //fqLens {1, 2, 4, 8, 16};
  {
    int node = 8;
    int fqLen = 2;
    Var fq, nd;

    Tuple a = lvl0_2(x, 2 * nd);
    Tuple b = lvl0_2(x, 2 * nd + 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl1_2(x, fq, nd) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl1_2.compute_root();
    lvl1_2.vectorize(nd, 8);
    // lvl1.trace_stores();
    // lvl1.realize({2, 8});
  }
  //lvl2
  {
    int node = 4;
    int fqLen = 4;
    Var fq, nd;

    Tuple a = lvl1_2(x, fq % fqLens(1), 2 * nd);
    Tuple b = lvl1_2(x, fq % fqLens(1), 2 * nd + 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    // Tuple mult = {b * e[0], b * e[1]};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl2_2(x, fq, nd) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl2_2.compute_root();
    // lvl2.trace_stores();
    // lvl2.realize({4, 4});
  }
  //lvl3
  {
    int node = 2;
    int fqLen = 8;
    // RDom fq(0, fqLen);
    // RDom nd(0, node);
    Var fq, nd;

    Tuple a = lvl2_2(x, fq % fqLens(2), 2 * nd);
    Tuple b = lvl2_2(x, fq % fqLens(2), 2 * nd + 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};


    lvl3_2(x, fq, nd) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl3_2.compute_root();
    // lvl3.trace_stores();
  }
  //lvl4
  {
    int node = 1;
    int fqLen = 16;
    // RDom fq(0, fqLen);
    Var fq;

    Tuple a = lvl3_2(x, fq % fqLens(3), 0);
    Tuple b = lvl3_2(x, fq % fqLens(3), 1);

    Tuple e = {cos(2.f * PI * fq / fqLen), sin(2.f * PI * fq / fqLen)};
    Tuple mult = {b[0] * e[0] - b[1] * e[1], b[0] * e[1] + b[1] * e[0]};

    lvl4_2(x, fq) = Tuple(a[0] + mult[0], a[1] + mult[1]);
    lvl4_2.compute_root();
    cout << endl << "print lvl4" << endl;
    lvl4_2.trace_stores();
  }
  out(x, y) = lvl4_2(x, y)[0] / 256.f;//(16.f * 16.f);
  
  return out;
}




void hft2D()
{
  Buffer<float> funB(16,16);
  for (int i = 0; i < 16; i++)
    for (int x = 0; x < 16; x++)
      funB(x, i) = float(i + x);
  
  Func fun;
  Var x, y;
  fun(x, y) = funB(x, y);

  Buffer<float> fft = fftAlg2D(fun).realize({16, 16});
  
  // cout << "ft: " << endl;
  // for (int i = 0; i < N; i++)
  // {
  //   cout << fft(0, i) << "   +    ";
  //   cout << fft(1, i) << "i" << endl;
  //   // cout << fft(i).real << "   +    ";
  //   // cout << fft(i).imag << "i" << endl;
  // }

  Func ftFun = fftAlg2D(fun);
  Func ftF;
  ftF(x, y) = ftFun(x, y);
  // ftF(x, y) = fft(x, y);
  // ftF(x) = Tuple(fft(0, x), fft(1, x));
  Buffer<float> ifft = ifftAlg2D(ftF).realize({16, 16});

  cout << "ift: " << endl;
  for (int i = 0; i < 16; i++)
    for (int x = 0; x < 16; x++)
    cout << ifft(x, i) << endl;
}


int main()
{
  hft2D();
  // ftN();
  
  
  
  
}