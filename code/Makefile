#Makefile 

#CC = icc
CC = g++

all: graph_matching gen_look_up

gen_look_up: gen_look_up.cpp
	$(CC) gen_look_up.cpp -o ./gen_look_up -pthread

graph_matching: main.cpp
	$(CC) main.cpp -o ./graph_matching

clean:
	rm -f ./graph_matching *.o ./gen_look_up
