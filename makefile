OS := $(shell uname)
ifeq ($(OS),Darwin)
        #CC      = /usr/local/opt/llvm/bin/clang++
        #CFLAGS  = -O3 -mavx -std=c++14 -w -march=native -I/usr/local/opt/llvm/include -fopenmp 
        #LDFLAGS = -L/usr/local/opt/llvm/lib
	CC      = g++
        CFLAGS  = -O3 -mavx -std=c++14 -w -march=native -I/usr/local/opt/libomp/include -fopenmp
        LDFLAGS = -L/usr/local/opt/libomp/lib
else
        CC      = g++
        CFLAGS  = -O3 -std=c++14 -w -fopenmp 
        LDFLAGS =
endif


SOURCES = containers/relation.cpp
OBJECTS = $(SOURCES:.cpp=.o)
	

all: main

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) main.cpp -o sj $(LDADD)

debug: $(OBJECTS)
	$(CC) -g $(CFLAGS) $(LDFLAGS) $(OBJECTS) main.cpp -o sj $(LDADD)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

.cc.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf containers/*.o
	rm -rf algorithms/*.o
	rm -rf partitioning/*.o
	rm -rf scheduling/*.o
	rm -rf sj