build:
	mpicc filtru.c -o filtru -lm

clean:
	rm -f filtru