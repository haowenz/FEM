c_source=sequence_batch.c index.c filter.c align.c input_queue.c output_queue.c map.c FEM_map.c FEM_index.c FEM.c
src_dir=src
objs_dir=objs
objs+=$(patsubst %.c,$(objs_dir)/%.o,$(c_source))

cxx=gcc
cxxflags=-g -Wall -O3 -march=native

ldflags=-lpthread -lm -lz
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
