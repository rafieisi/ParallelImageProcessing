#include "pgm.h"
#include <stdio.h>
#include <stdlib.h>

int main(){
pgm_image image1,image2;
    create_random_pgm_image(&image1, 100000000, 1);
    save_pgm_to_file("1row.pgm", &image1);

    create_random_pgm_image(&image2, 1, 100000000);
    save_pgm_to_file("1column.pgm", &image2);
}