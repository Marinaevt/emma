#pragma once

// TODO: написать свои реализации представленных функций, ибо использовать ворованый код Iowa-Hills P51 root finder не комильфо
// http://iowahills.com/P51RootFinder.html

long double BuildPolyFromRoots(int N, long double *RealRoot, long double *ImagRoot, long double *PolyCoeff);
void SortRoots(int N, long double *OrigRealRoot, long double *OrigImagRoot, long double *RealRoot, long double *ImagRoot);

int QuadCubicRoots(long double *Coeff, int N, long double *RealRoot, long double *ImagRoot);
void QuadRoots(long double *P, long double *RealPart, long double *ImagPart);
void CubicRoots(long double *P, long double *RealPart, long double *ImagPart);
void BiQuadRoots(long double *P, long double *RealPart, long double *ImagPart);
void ReversePoly(long double *P, int N);
void InvertRoots(int N, long double *RealRoot, long double *ImagRoot);

inline void ShowMessage(const char* str) {}