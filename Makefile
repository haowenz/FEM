c_source=utils.c readtools.c reftools.c indextools.c filter.c verifier.c outputer.c mapper.c FEM-align.c FEM-index.c FEM.c
src_dir=src
objs_dir=objs
objs+=$(patsubst %.c,$(objs_dir)/%.o,$(c_source))

cxx=gcc
cxxflags=-Wall -O3 -march=native

ldflags=-lpthread -lm
exec=FEM

all: dir $(exec) 
	
dir:
	-mkdir -p $(objs_dir)

$(exec): $(objs)
	$(cxx) $(cxxflags) $(objs) -o $(exec) $(ldflags)
	
$(objs_dir)/%.o: $(src_dir)/%.c
	$(cxx) $(cxxflags) -c $< -o $@

.PHONY: clean
clean:
	-rm -r $(exec) $(objs_dir)
