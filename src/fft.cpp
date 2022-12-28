#include "fft.hpp"

/*
https://algoteka.com/samples/45/polynomial-multiplication-c-plus-plus-o%2528n-log-n%2529-fast-fourier-transform-based-solution
*/

std::vector<std::complex<double>> fast_fourier_transform(std::vector<std::complex<double>> x, bool inverse) {
    std::vector<std::complex<double>> w(x.size(), 0.0);

    w[0] = 1.0;
    for(auto pow_2 = 1; pow_2 < x.size(); pow_2 *= 2) {
        w[pow_2] = std::polar(1.0, 2*M_PI*pow_2/x.size() * (inverse ? 1 : -1));
    }
    for(auto i=3, last=2; i < x.size(); ++i) {
        if(w[i] == 0.0) {
            w[i] = w[last]*w[i-last];
        }
        else {
            last = i;
        }
    }

    for(auto block_size = x.size(); block_size > 1; block_size /= 2) {
        std::vector<std::complex<double>> new_x(x.size());

        for(auto start = 0; start < x.size(); start += block_size) {
            for(int i = 0; i < block_size; ++i) {
                new_x[start+block_size/2 * (i%2) + i/2] = x[start + i];
            }
        }
        x = new_x;
    }

    for(auto block_size = 2; block_size <= x.size(); block_size *= 2) {
        std::vector<std::complex<double>> new_x(x.size());

        int w_base_i = x.size() / block_size;

        for(auto start = 0; start < x.size(); start += block_size) {
            for(auto i = 0; i < block_size/2; ++i) {
                new_x[start + i] = x[start+i] + w[w_base_i*i] * x[start + block_size/2 + i];
                new_x[start+block_size/2+i] = x[start+i] - w[w_base_i*i] * x[start + block_size/2 + i];
            }
        }
        x = new_x;
    }

    return x;
}

std::vector<std::complex<double>> fft(std::vector<std::complex<double>> x) {
    return fast_fourier_transform(x, false);
}

std::vector<std::complex<double>> ifft(std::vector<std::complex<double>> x) {
    return fast_fourier_transform(x, true);
}