AM_CFLAGS = -O2
AM_CPPFLAGS = -DDATA_PATH=\"$(pkgdatadir)\"

bin_PROGRAMS = mia ma ccheck

mia_SOURCES = mia.c mia.h params.h types.h pssm.c pssm.h fsdb.h fsdb.c kmer.c kmer.h mia_main.c map_align.c map_align.h io.h io.c map_alignment.h map_alignment.c

ma_SOURCES = params.h types.h map_alignment.h map_alignment.c map_assembler.c io.h io.c map_align.h map_align.c

ccheck_SOURCES = ccheck.cc myers_align.c fsdb.c io.c kmer.c map_align.c map_alignment.c mia.c pssm.c mt311.c \
		 map_align.h params.h types.h io.h map_alignment.h config.h mia.h fsdb.h pssm.h kmer.h myers_align.h
