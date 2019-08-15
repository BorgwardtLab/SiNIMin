all: base
base:
	g++ -O5 -std=c++11 -march=native -I/usr/local/include/eigen3 \
		-o sinimin ../C/main/edge_epistasis.cpp -fopenmp

	g++ -O5 -std=c++11 -march=native -I/usr/local/include/eigen3 \
		-o sinimin_fwer ../C/main/sinimin_fwer_estimation.cpp -fopenmp


