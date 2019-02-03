
# Note to CS149 students: change the REF_LIB define depending on the
# platform you are running on.  Different precompiled reference
# libraries are needed for testing on myth and GCP

# REF_LIB := ref_pr_myth.a
REF_LIB := ref_pr_cloud.a


all: default grade

default: page_rank.cpp main.cpp
	g++ -I../ -std=c++11 -fopenmp -g -O3 -o pr main.cpp page_rank.cpp ../common/graph.cpp ${REF_LIB}
grade: page_rank.cpp grade.cpp
	g++ -I../ -std=c++11 -fopenmp -g -O3 -o pr_grader grade.cpp page_rank.cpp ../common/graph.cpp ${REF_LIB}
clean:
	rm -rf pr pr_grader *~ *.*~
