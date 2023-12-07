SRC = main.cc run.cc
all: $(SRC)
	mpicxx $(SRC) -o matmul 
.PHONY: clean
clean:
	rm -rf matmul