all: base
base:
	g++-6 -O5 -std=c++11 -march=native -I/usr/local/include/eigen3 \
		-o ../compiled/sinimin ./main/sinimin.cpp -fopenmp

	g++-6 -O5 -std=c++11 -march=native -I/usr/local/include/eigen3 \
		-o ../compiled/sinimin_fwer ./main/sinimin_fwer_estimation.cpp -fopenmp


