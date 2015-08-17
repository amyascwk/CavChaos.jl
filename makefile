all:
	gcc -c -Wall -Werror -fPIC -o bin/ray.o src/ray.c
	gcc -c -Wall -Werror -fPIC -o bin/mode.o src/mode.c
	gcc -shared -o lib/libcavchaos.so bin/ray.o bin/mode.o -lgsl

