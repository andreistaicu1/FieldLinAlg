#ifndef FFT_H
#define FFT_H
#include <complex>
#include <vector>
#include <cmath>

extern std::vector<std::complex<double>> fft(std::vector<std::complex<double>>);
extern std::vector<std::complex<double>> ifft(std::vector<std::complex<double>>);

#endif