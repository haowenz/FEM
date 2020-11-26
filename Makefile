c_source=sequence_batch.c index.c filter.c align.c input_queue.c output_queue.c map.c FEM_map.c FEM_index.c FEM.c
src_dir=src
objs_dir=objs
objs+=$(patsubst %.c,$(objs_dir)/%.o,$(c_source))

fem_dir = $(shell pwd)
htslib_dir ?= ${fem_dir}/extern/htslib/build
htslib_include_dir ?= ${htslib_dir}/include
htslib_lib_dir ?= ${htslib_dir}/lib
htslib_lib ?= hts

cxx=gcc
cxxflags=-Wall -O3 -march=native -I${htslib_include_dir}

ldflags=-L${htslib_lib_dir} -l${htslib_lib} -lpthread -lm -lz
exec=FEM

.PHONY: all
all: htslib check_htslib mk_obj_dir $(exec) 

.PHONY: htslib
htslib:
	git submodule update --init --recursive
	cd extern/htslib;\
	mkdir build;\
	autoheader;\
	autoconf;\
	./configure CC="${cxx}" --disable-bz2 --disable-lzma --prefix="${htslib_dir}";\
	make -j;\
	make install

.PHONY: check_htslib
check_htslib:
	@[ -f "${htslib_include_dir}/htslib/sam.h" ] || { echo "htslib headers not found" >&2; exit 1; }
	@[ -f "${htslib_lib_dir}/lib${htslib_lib}.so" ] || [ -f "${htslib_lib_dir}/lib${htslib_lib}.a" ] || { echo "htslib library not found" >&2; exit 1; }

.PHONY: mk_obj_dir
mk_obj_dir:
	mkdir -p $(objs_dir)

$(exec): $(objs)
	$(cxx) $(cxxflags) $(objs) -o $(exec) $(ldflags)
	
$(objs_dir)/%.o: $(src_dir)/%.c
	$(cxx) $(cxxflags) -c $< -o $@

.PHONY: clean
clean:
	cd "extern/htslib" && make clean
	rm -rf "${htslib_dir}"
	rm -rf $(exec) $(objs_dir)
