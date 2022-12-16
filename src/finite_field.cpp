#include <vector>
#include <memory>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <iostream>

#include "min_polynomial.cpp"

#include <tuple>

class Polynomial {
public:
    std::vector<int> p;
    int actual_degree;

    Polynomial(std::vector<int> p) : p(p) {
        int last = 0;
        for(auto i = 0; i < this->p.size(); ++i) {
            if(this->p[i] != 0) {
                last = i;
            }
        }
        this->actual_degree = last;
    };

    bool operator==(const Polynomial& other) {

        if (this->actual_degree != other.actual_degree){
            return false;
        }

        // Gets the smallest length polynomial
        int min = this->p.size();
        if (other.p.size() < min){
            min = other.p.size();
        }

        // Checks each element to see if they're equal
        for (auto i = 0; i < min; i++) {
            if (this->p[i] != other.p[i]) return false;
        }

        return true;
    }

    Polynomial operator+(const Polynomial& other) {
        auto max = std::max(this->p.size(), other.p.size());

        std::vector<int> output(max,0);
        for(auto i = 0; i < max; ++i) {
            int first = 0;
            if(this->p.size() > i) {
                first = this->p[i];
            }
            
            int second = 0;
            if(other.p.size() > i) {
                second = other.p[i];
            }

            output[i] = first+second;
        }

        return Polynomial(output);
    };

    Polynomial operator-(const Polynomial& other) {
        auto max = std::max(this->p.size(), other.p.size());

        std::vector<int> output(max,0);
        for(auto i = 0; i < max; ++i) {
            int first = 0;
            if(this->p.size() > i) {
                first = this->p[i];
            }
            
            int second = 0;
            if(other.p.size() > i) {
                second = other.p[i];
            }

            output[i] = first-second;
        }

        return Polynomial(output);
    };

    Polynomial operator*(const Polynomial& other){
        
        int new_degree = this->actual_degree + other.actual_degree + 1;
        std::vector<int> output(new_degree, 0);

        for(auto n = 0; n < new_degree; n++){
            for(auto i = 0; i <= n; i++){
                if (i <= this->actual_degree && n - i <= other.actual_degree){
                    output[n] += this->p[i] * other.p[n-i];
                }
            }
        }
        
        return Polynomial(output);
    }
    
    Polynomial trim_end() {
        std::vector<int> output;
        for(auto i = 0; i < actual_degree+1; ++i) {
            output.push_back(this->p[i]);
        }

        return Polynomial(output);
    }

    Polynomial pad_to_length(int length) {
        std::vector<int> output;
        for(auto i = 0; i < length; ++i) {
            if(i < this->p.size()) {
                output.push_back(this->p[i]);
            }
            else {
                output.push_back(0);
            }
        }
        return Polynomial(output);
    }

    int size() {
        return this->p.size();
    }

    int lead_coeff() {
        return this->p[this->p.size()-1];
    }

    void print() {
        for(auto elem : this->p) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    static Polynomial zero_polynomial(int length) {
        std::vector<int> output(length, 0);
        return Polynomial(output);
    }

    static std::tuple<Polynomial, Polynomial> polynomial_div(Polynomial p1, Polynomial p2) {
        int max_size = std::max(p1.size(), p2.size());
        Polynomial zero_poly = zero_polynomial(max_size);

        if(p2.actual_degree > p1.actual_degree) {
            return std::make_pair(zero_poly, p1.pad_to_length(max_size));
        }

        Polynomial trim_p1 = p1.trim_end();
        Polynomial trim_p2 = p2.trim_end();

        int coefficient = trim_p1.lead_coeff()/trim_p2.lead_coeff();
        int degree = trim_p1.size()-1 - (trim_p2.size()-1);
        std::vector<int> on_top_vec(degree+1, 0);
        on_top_vec[degree] = coefficient;
        Polynomial on_top = Polynomial(on_top_vec);

        Polynomial remainder = (trim_p1 - (on_top*trim_p2)).trim_end();
        
        std::tuple<Polynomial, Polynomial> next_output = polynomial_div(remainder, trim_p2);
        return std::make_pair((on_top + std::get<0>(next_output)).pad_to_length(max_size), std::get<1>(next_output).pad_to_length(max_size));
    }

    static Polynomial euclidean_algorithm(Polynomial p1, Polynomial p2, std::tuple<Polynomial, Polynomial> compute_p1, 
        std::tuple<Polynomial, Polynomial> compute_p2){

            auto output = polynomial_div(p1, p2);
            Polynomial quotient = std::get<0>(output);
            Polynomial remainder = std::get<1>(output);

            // This function should only be input relatively prime things
            assert(!(remainder == zero_polynomial(1)));

            std::vector<int> one_vec{1};
            Polynomial one = Polynomial(one_vec);

            // Checks if the remainder is one
            if (remainder == one){
                return std::get<1>(compute_p1) - (std::get<1>(compute_p2) * quotient);
            }

            Polynomial compute_rem_a = std::get<0>(compute_p1) - (std::get<0>(compute_p2) * quotient);
            Polynomial compute_rem_b = std::get<1>(compute_p1) - (std::get<1>(compute_p2) * quotient);
            return euclidean_algorithm(p2, remainder, compute_p2, std::make_pair(compute_rem_a, compute_rem_b));
    }
};

int main() {
    std::vector<int> p2_vec{1,1,1};
    std::vector<int> p1_vec{2,0,0,2,5,6};
    Polynomial p2 = Polynomial(p2_vec);
    Polynomial p1 = Polynomial(p1_vec);

    auto output = Polynomial::polynomial_div(p1,p2);

    std::cout << "Quotient: ";
    for(int i = 0; i < std::get<0>(output).size(); ++i) {
        std::cout << std::get<0>(output).p[i] << " ";
    } 
    std::cout << std::endl;

    std::cout << "Remainder: ";
    for(int i = 0; i < std::get<1>(output).size(); ++i) {
        std::cout << std::get<1>(output).p[i] << " ";
    } 
    std::cout << std::endl;
    
    return 0;
}


class finite_field {
public:
    int p_char, degree;
    std::vector<int> min_polynomial;
    
    finite_field(int p_char, int degree) : p_char(p_char), degree(degree) {
        std::vector<int> temp(10, 0);
        this->min_polynomial = temp;
    };

    bool operator== (finite_field const &other) {
        return this->p_char == other.p_char && this->degree == other.degree;
    }    
    
};

class ff_element {
public:
    std::vector<int> polynomial;
    finite_field ff;

    ff_element(std::vector<int> polynomial, finite_field ff) : polynomial(polynomial), ff(ff) {};
    ff_element(finite_field ff) : ff(ff) {
        std::vector<int> zeroes(ff.degree, 0);
        this->polynomial = zeroes;
    };

    ff_element operator+ (ff_element const &other) {
        assert(this->ff == other.ff);

        ff_element output = ff_element(this->ff);
        for (auto i = 0; i < this->ff.degree; i++){
            output.polynomial[i] = (this->polynomial[i]+other.polynomial[i]) % this->ff.p_char;
        }
        
        return output;
    };

    ff_element operator- (ff_element const &other) {
        assert(this->ff == other.ff);
        
        ff_element output = ff_element(this->ff);
        for (auto i = 0; i < this->ff.degree; i++){
            output.polynomial[i] = (this->polynomial[i] - other.polynomial[i]) % this->ff.p_char;
        }

        return output;
    };

    ff_element operator* (ff_element const &other) {
        assert(this->ff == other.ff);
        
        ff_element output = ff_element(this->ff);
    };

    ff_element operator/ (ff_element const &other) {
        assert(this->ff == other.ff);
        
        ff_element output = ff_element(this->ff);
        
    };
};

