TARGET= aos.out aos_pair.out aos_intrin.out soa.out soa_pair.out soa_intrin.out aos_intrin_mat_transpose.out aos_sorted_z.out 

CC=g++
CPPFLAGS=-O3 -std=c++11 -march=native

all: $(TARGET)

aos.out: force_aos.cpp
	$(CC) $(CPPFLAGS) $< -o $@

aos_pair.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DPAIR $< -o $@

aos_intrin.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DINTRIN $< -o $@

aos_sorted_z.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DSORTED_Z $< -o $@

aos_sorted_z_intrin.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DSORTED_Z_INTRIN $< -o $@

aos_intrin_mat_transpose.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DMAT_TRANSPOSE $< -o $@

soa.out: force_soa.cpp
	$(CC) $(CPPFLAGS) $< -o $@

soa_pair.out: force_soa.cpp
	$(CC) $(CPPFLAGS) -DPAIR $< -o $@

soa_intrin.out: force_soa.cpp
	$(CC) $(CPPFLAGS) -DINTRIN $< -o $@

clean:
	rm -f $(TARGET)

test: aos_pair.out aos_intrin.out soa_pair.out soa_intrin.out aos_intrin_mat_transpose.out
	./aos_pair.out > aos_pair.dat
	./aos_intrin.out > aos_intrin.dat
	./aos_intrin_mat_transpose.out > aos_intrin_mat_transpose.dat
	diff aos_pair.dat aos_intrin.dat
	diff aos_intrin.dat aos_intrin_mat_transpose.dat
	./soa_pair.out > soa_pair.dat
	./soa_intrin.out > soa_intrin.dat
	diff soa_pair.dat soa_intrin.dat

-include makefile.opt
