c_source=utils.c readtools.c reftools.c indextools.c filter.c verifier.c outputer.c mapper.c FEM-align.c FEM-index.c FEM.c
src_dir=src
objs_dir=objs
objs+= $(patsubst %.c, $(objs_dir)/%.o, $(c_source))

cxx=g++
cxxflags =-std=c++11 -Wall -O3 -funroll-all-loops -fopenmp -march=native -lpthread

ldflags= 
exec=FEM

all: dir $(exec) 
	
dir:
	-mkdir -p $(objs_dir)

$(exec): $(objs)
	$(cxx) $(cxxflags) $(objs) -o $(exec)
	
clean:
	-rm -rf $(objs)

$(objs_dir)/%.o: $(src_dir)/%.c
	$(cxx) $(cxxflags) -o $@ -c $< 
