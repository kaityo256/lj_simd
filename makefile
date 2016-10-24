TARGET= a.out pair.out intrin.out soa.out
all: $(TARGET)

a.out: test.cpp
	icpc -O3 -xHOST -std=c++11 $< -o $@

soa.out: test_soa.cpp
	icpc -O3 -xHOST -std=c++11 $< -o $@

pair.out: test.cpp
	icpc -O3 -xHOST -std=c++11 -DPAIR $< -o $@

intrin.out: test.cpp
	icpc -O3 -xHOST -std=c++11 -DINTRIN $< -o $@

clean:
	rm -f $(TARGET)

test:
	./pair.out > pair.dat
	./intrin.out > intrin.dat
	diff pair.dat intrin.dat

