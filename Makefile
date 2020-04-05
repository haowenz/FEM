c_source=sequence_batch.c index.c filter.c align.c input_queue.c output_queue.c map.c FEM_map.c FEM_index.c FEM.c
src_dir=src
objs_dir=objs
objs+=$(patsubst %.c,$(objs_dir)/%.o,$(c_source))

cxx=gcc
cxxflags=-g -Wall -O3 -march=native

ldflags=-lpthread -lm -lz
exec=FEM

.PHONY: all
all: htslib mk_obj_dir $(exec) 

.PHONY: htslib
htslib:
	git submodule update --init --recursive
	cd extern/htslib;\
	autoheader;\
	autoconf;\
	./configure --disable-bz2 --disable-lzma;\
	make;\
	make prefix=. install

.PHONY: mk_obj_dir
mk_obj_dir:
	mkdir -p $(objs_dir)

$(exec): $(objs)
	$(cxx) $(cxxflags) $(objs) -o $(exec) $(ldflags)
	
$(objs_dir)/%.o: $(src_dir)/%.c
	$(cxx) $(cxxflags) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(exec) $(objs_dir)
