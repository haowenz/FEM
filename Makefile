c_source=utils.c readtools.c reftools.c indextools.c filter.c verifier.c outputer.c mapper.c FEM-align.c FEM.c
cpp_source=indexer.cc FEM-index.cc
src_dir=src
objs_dir=objs
objs+= $(patsubst %.c, $(objs_dir)/%.o, $(c_source))
objs+= $(patsubst %.cc, $(objs_dir)/%.o, $(cpp_source))

cxx=g++
cxxflags =-std=c++11 -O3 -funroll-all-loops -fopenmp -march=native -lpthread

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

$(objs_dir)/%.o: $(src_dir)/%.cc
	$(cxx) $(cxxflags) -o $@ -c $<
