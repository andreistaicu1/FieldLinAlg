#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "finite_field.cpp"
#include "polynomial.cpp"

int normal_main(){
    int p = 23;
    std::vector<int> values = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                            14, 15, 16, 17, 18, 19, 20, 21, 22};
    
    finite_field ff = finite_field(p, 1);

    for (int val : values){
        std::cout << ff.inverse_mod_p(val) << " ";
    }
    std::cout << std::endl;

    return 0;
}

int main(int argc, char** argv) {
    doctest::Context context;
    context.applyCommandLine(argc, argv);
    int res = context.run();
    if(context.shouldExit()) {
        return res;
    }

    int return_code = normal_main();

    return return_code+res;
}