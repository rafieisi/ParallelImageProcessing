CC=gcc
GCC_OPT = -O2 -Wall -Werror

%.o: %.c
	$(CC) -c -o $@ $< $(GCC_OPT)

main: very_big_sample.o very_tall_sample.o main.c pgm.c filters.c
	$(CC) $(GCC_OPT) main.c pgm.c filters.c very_big_sample.o very_tall_sample.o -o main.out -lpthread

#TODO: run your script that runs the experiments, collects all necessary data, and produces the graphs
run: 
	 gcc imagecreator.c pgm.c -o imagecreator
	 ./imagecreator
	 python3 perfs_student.py


pgm_creator:
	$(CC) $(GCC_OPT) pgm_creator.c pgm.c -o pgm_creator.out
	
clean:
	rm *.o *.out very* *.pgm imagecreator perfb* *.png *.pickle rm -r *pycache*
