CC=g++
CPPFLAGS=-O3 -std=c++11 -march=native

-include makefile.opt

ifdef AVX512
DAVX512=-DAVX512
endif

AOS_BINS= aos.out aos_pair.out aos_sorted.out
AOS_BINS +=aos_avx2.out
AOS_BINS +=aos_avx2_swp.out
ifdef AVX512
AOS_BINS +=aos_avx512.out
AOS_BINS +=aos_avx512_loopopt.out
AOS_BINS +=aos_avx512_loopopt_swp.out
AOS_BINS +=aos_avx512_gatheronly.out
AOS_BINS +=aos_avx512_transpose.out
endif

SOA_BINS =soa.out soa_pair.out soa_sorted.out
SOA_BINS +=soa_avx2.out
ifdef AVX512
SOA_BINS +=soa_avx512.out
SOA_BINS +=soa_avx512_loopopt.out
SOA_BINS +=soa_avx512_loopopt_swp.out
endif

TARGET = $(AOS_BINS) $(SOA_BINS)


all: $(TARGET)

aos.out: force_aos.cpp
	$(CC) $(CPPFLAGS) $(DAVX512) $< -o $@

aos_pair.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DPAIR $< -o $@

aos_sorted.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DSORTED $< -o $@

aos_avx2.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DAVX2 $< -o $@

aos_avx2_swp.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DAVX2_SWP $< -o $@

aos_sorted_z.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DSORTED_Z $< -o $@

aos_sorted_z_avx2.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DSORTED_Z_AVX2 $< -o $@

aos_avx512.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DAVX512 -DAVX512_SIMPLE $< -o $@

aos_avx512_loopopt.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DAVX512 -DAVX512_LOOPOPT $< -o $@

aos_avx512_loopopt_swp.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DAVX512 -DAVX512_LOOPOPT_SWP $< -o $@

aos_avx512_gatheronly.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DAVX512 -DAVX512_GATHERONLY $< -o $@

aos_avx512_transpose.out: force_aos.cpp
	$(CC) $(CPPFLAGS) -DAVX512 -DAVX512_TRANSPOSE $< -o $@

soa.out: force_soa.cpp
	$(CC) $(CPPFLAGS) $(DAVX512) $< -o $@

soa_pair.out: force_soa.cpp
	$(CC) $(CPPFLAGS) -DPAIR $< -o $@

soa_sorted.out: force_soa.cpp
	$(CC) $(CPPFLAGS) -DSORTED $< -o $@

soa_avx2.out: force_soa.cpp
	$(CC) $(CPPFLAGS) -DAVX2 $< -o $@

soa_avx512.out: force_soa.cpp
	$(CC) $(CPPFLAGS) -DAVX512 -DAVX512_SIMPLE $< -o $@

soa_avx512_loopopt.out: force_soa.cpp
	$(CC) $(CPPFLAGS) -DAVX512 -DAVX512_LOOPOPT $< -o $@

soa_avx512_loopopt_swp.out: force_soa.cpp
	$(CC) $(CPPFLAGS) -DAVX512 -DAVX512_LOOPOPT_SWP $< -o $@


clean:
	rm -f $(TARGET)

test: $(TARGET)
	./aos_pair.out > aos_pair.dat
	./aos_avx2.out > aos_avx2.dat
	diff aos_pair.dat aos_avx2.dat
	./soa_pair.out > soa_pair.dat
	./soa_avx2.out > soa_avx2.dat
	diff soa_pair.dat soa_avx2.dat

REF=aos_avx2.out
SUB=aos_avx2_swp.out

test2: $(REF) $(SUB)
	./$(REF) > orig.dat
	./$(SUB) > test.dat
	diff orig.dat test.dat

