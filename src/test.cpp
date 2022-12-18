#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "finite_field.cpp"
#include "polynomial.cpp"

int main(){
    int p = 23;
    std::vector<int> values = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                            14, 15, 16, 17, 18, 19, 20, 21, 22};
    
    finite_field ff = finite_field(p, 1);

    for (int val : values){
        std::cout << ff.inverse_mod_p(val) << " ";
    }
    std::cout << std::endl;
}