TARGET= aos.out aos_pair.out aos_intrin.out soa.out soa_pair.out soa_intrin.out

all: $(TARGET)

aos.out: force_aos.cpp
	icpc -O3 -xHOST -std=c++11 $< -o $@

aos_pair.out: force_aos.cpp
	icpc -O3 -xHOST -std=c++11 -DPAIR $< -o $@

aos_intrin.out: force_aos.cpp
	icpc -O3 -xHOST -std=c++11 -DINTRIN $< -o $@

soa.out: force_soa.cpp
	icpc -O3 -xHOST -std=c++11 $< -o $@

soa_pair.out: force_soa.cpp
	icpc -O3 -xHOST -std=c++11 -DPAIR $< -o $@

soa_intrin.out: force_soa.cpp
	icpc -O3 -xHOST -std=c++11 -DINTRIN $< -o $@

clean:
	rm -f $(TARGET)

test: aos_pair.out aos_intrin.out soa_pair.out soa_intrin.out
	./aos_pair.out > aos_pair.dat
	./aos_intrin.out > aos_intrin.dat
	diff aos_pair.dat aos_intrin.dat
	./soa_pair.out > soa_pair.dat
	./soa_intrin.out > soa_intrin.dat
	diff soa_pair.dat soa_intrin.dat

