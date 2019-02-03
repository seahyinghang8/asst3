
# Note to CS149 students: change the REF_LIB define depending on the
# platform you are running on.  Different precompiled reference
# libraries are needed for testing on myth and GCP

REF_LIB :=  ref_bfs_myth.o
#REF_LIB :=  ref_bfs_cloud.o


all: default grade

default: main.cpp bfs.cpp
	g++ -I../ -std=c++11 -fopenmp -O3 -g -o bfs main.cpp bfs.cpp ../common/graph.cpp ${REF_LIB}
grade: grade.cpp bfs.cpp
	g++ -I../ -std=c++11 -fopenmp -O3 -g -o bfs_grader grade.cpp bfs.cpp ../common/graph.cpp ${REF_LIB}
clean:
	rm -rf bfs_grader bfs  *~ *.*~
