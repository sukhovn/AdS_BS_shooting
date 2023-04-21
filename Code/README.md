The modes section of this code requires the following libraries:
Eigen 3
FFTW 3

To install Eigen 3 download the latest release from https://eigen.tuxfamily.org/,
copy the folder Eigen from the archive to /usr/local/include/eigen3/

To install FFTW 3 download the latest release from https://www.fftw.org/download.html,
extract the archive contents and run:
	./configure CFLAGS="-fPIC"
	make
	make install