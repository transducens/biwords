# author: Felipe Sánchez Martínez

bin_PROGRAMS = symmetrize-alignments gen-text-to-compress \
               stats-alignments restore-text filter-alignments \
               restore-text-multiple-offsets \
               count-unaligned-words

bin_SCRIPTS = giza++A3-to-txt switch-alignments get-one-side

#bitext_compressiondir=$(prefix)/bitext-compression-$(GENERIC_API_VERSION)
#bitext_compression_DATA = apertium-transfer-at.dtd apertium-posttransfer.dtd

EXTRA_DIST = giza++A3-to-txt.awk switch-alignments.awk get-one-side.pl

symmetrize_alignments_SOURCES = symmetrize-alignments.cc \
       alignment.h alignment.cc utils.h utils.cc zfstream.h zfstream.cc 
symmetrize_alignments_LDADD = $(Z_LIBS)

stats_alignments_SOURCES = stats-alignments.cc \
       alignment.h alignment.cc utils.h utils.cc zfstream.h zfstream.cc 
stats_alignments_LDADD = $(Z_LIBS)

gen_text_to_compress_SOURCES = gen-text-to-compress.cc \
       alignment.h alignment.cc utils.h utils.cc zfstream.h zfstream.cc 
gen_text_to_compress_LDADD = $(Z_LIBS)

restore_text_multiple_offsets_SOURCES = restore-text-multiple-offsets.cc utils.h utils.cc alignment.h alignment.cc

restore_text_SOURCES = restore-text.cc utils.h utils.cc alignment.h alignment.cc

filter_alignments_SOURCES = filter-alignments.cc \
       alignment.h alignment.cc zfstream.h zfstream.cc utils.h utils.cc
filter_alignments_LDADD = $(Z_LIBS)

count_unaligned_words_SOURCES = count-unaligned-words.cc \
       alignment.h alignment.cc zfstream.h zfstream.cc utils.h utils.cc
count_unaligned_words_LDADD = $(Z_LIBS)

INCLUDES = $(Z_CFLAGS)

CLEANFILES = *~ $(bin_SCRIPTS)

giza++A3-to-txt:  giza++A3-to-txt.awk Makefile.am
	@echo "Creating $@ script"
	@echo "#!$(AWK) -f" > $@
	@cat giza++A3-to-txt.awk >> $@
	@chmod a+x $@

switch-alignments: switch-alignments.awk Makefile.am
	@echo "Creating $@ script"
	@echo "#!$(AWK) -f" > $@
	@cat switch-alignments.awk >> $@
	@chmod a+x $@
	
get-one-side: get-one-side.pl Makefile.am
	@echo "Creating $@ script"
	@echo "#!$(PERL) -w" > $@  
	@cat get-one-side.pl >> $@
	@chmod a+x $@
