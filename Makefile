EXPORTDIR = ../C++tree.export
WEBDIR = /usr2/makino/WWW/softwares/C++tree
PGPLOT_DIR = /usr/local/pgplot
#PGPLOT_DIR = /data6/makino/PDS/pgplot_linux/pgplot
PGINCLUDE =  -I$(PGPLOT_DIR)
CFLAGS = $(FLOAT_OPTIONS) -DREAL=double -O4 -ffast-math -funroll-loops #-g  # -p

OPTFLAGS = -march=skylake-avx512 -fopenmp -ftree-vectorize -fopt-info-vec
CLANGOPTFLAGS =  -fopenmp -ftree-vectorize
#OPTFLAGS = -march=skylake
#
# Add -DNOGRAPHICS for installation without PGPLOT
#CPPFLAGS =   -DTREE $(PGINCLUDE) -p $(FLOAT_OPTIONS) -O4   -ffast-math -funroll-loops
CPPFLAGS = $(OPTFLAGS) -DNOGRAPHICS -DTREE $(PGINCLUDE)  $(FLOAT_OPTIONS) -O4 -ffast-math -funroll-loops  #   -p
CLANGFLAGS = $(CLANGOPTFLAGS) -DNOGRAPHICS -DTREE $(PGINCLUDE)  $(FLOAT_OPTIONS) -O4 -ffast-math -funroll-loops  #  -pg -p
#CLANGFLAGS += -Rpass=.*
CLANGFLAGS += -Rpass=sve-loop-vectorize -Rpass-analysis=inline -Rpass-analysis=loop-vectorize
#PGLIB     =  -L$(PGPLOT_DIR) -g3 -O4 -lcpgplot -lpgplot -lX11 -lUfor -lfor -lFutil -lm -lprof1 -lpdf
PGLIB      = 
#PGLIB     =  -L$(PGPLOT_DIR) -lcpgplot -lpgplot -lm  -L/usr/X11R6/lib -lX11 -lg2c
#G4LIB =  -L/usr2/makino/src/pci-harp3/export/library -lg4a -lrt -lm
#G4LIB =  -L/usr6/kawai/pub/lib -lg4a -lrt -lm
G4LIB =  -L/usr2/makino/disk3src/harplibs -lg4a -lrt -lm
G6LIB =  -L/usr2/makino/disk3src/harplibs/x86_64 -lg6lx 
G6DIR = ../g6hib-x86

H3LIB = -L/usr2/makino/src/harp3board -lharp3  -lrt -lm 

#LIBOBJS = pgetopt.o gravity.o BHtree.o second.o
LIBOBJS = pgetopt.o gravity.o second.o
LIBOBJSE = pgetopt.o  BHtree.o second.o
LIBOBJSG4 = pgetopt.o  gravity_g4.o BHtree_g4.o second.o
LIBOBJSG6 = pgetopt.o  gravity_g6.o BHtree_g6.o second.o

EXPORT_FILES = BHtree.h nbody_particle.h vector.h \
               BHtree.C	gravity.C nbody.C pgetopt.C \
               nbtest.C\
	       second.c \
	       Makefile nbody_guide.tex\
               README COPYRIGHT \
               samplein sampleparm samplelog \
	       utils.hpp

ADDITIONAL_FILES = hom32.stoa hom64.stoa hom32k.stoa hom1M.stoa
A64FXFFLAGS        =  -O3 -cpp -g -pg -p -march=armv8-a+sve+simd # -fopenmp
#A64FLAGS      =  -Nfjprof,rt_tune_loop,rt_tune_func -Kocl,parallel,simd=2,optmsg=2,fast,loop_fission,swp_strong -O3\
#                -DUNIX -DBSD  -DFUJITSU
#A64FLAGS      =  -Nfjprof   -Khpctag,fast,ocl,simd=2,optmsg=2,loop_fission,swp_strong  -O3  -DUNIX -DBSD  -DFUJITSU
#A64FLAGS2      = -Nlst=t  -Kocl,simd=2,optmsg=2,fast -O3\
#                -DUNIX -DBSD  -DFUJITSU
A64FLAGS      =   -DNOGRAPHICS -DTREE -Ksimd=2   -Kfast -O3\
                -DUNIX -DBSD  -DFUJITSU
A64CCOMPILER   = FCC

FILESTOSENDTOFX700 = $(EXPORT_FILES) $(ADDITIONAL_FILES)
OTHERSRCS =  pgetopt.C gravity.C
OTHEROBJS =  second.o

pgetopt.cc: pgetopt.C
	cp -p  pgetopt.C pgetopt.cc
gravity.cc: gravity.C
	cp -p  gravity.C  gravity.cc
BHtree.cc: BHtree.C
	cp -p  BHtree.C BHtree.cc
nbtest.cc: nbtest.C
	cp -p  nbtest.C nbtest.cc




fx700: sendtofx700
sendtofx700: $(FILESTOSENDTOFX700)
	rsync -avp -e "ssh -l jmakino" $(FILESTOSENDTOFX700) login.cloud.r-ccs.riken.jp:src/C++tree

export: 
	rsync -avp $(EXPORT_FILES) $(EXPORTDIR)
	cd $(EXPORTDIR); git commit -a ; git push
export.tar.gz: $(EXPORT_FILES)
	tar cvf export.tar  $(EXPORT_FILES)
	gzip export.tar 
gravity.o : gravity.C  nbody_particle.h BHtree.h BHtree.C
	g++  $(CPPFLAGS)  -c gravity.C
gravity_g4.o : gravity.C  nbody_particle.h
	g++  -DHARP3 $(CPPFLAGS)  -c gravity.C
	mv gravity.o gravity_g4.o
gravity_g6.o : gravity.C  nbody_particle.h
	g++  -DGRAPE-6 $(CPPFLAGS)  -c gravity.C
	mv gravity.o gravity_g6.o



nbody.o : nbody.C  nbody_particle.h
	g++  -c    $(CPPFLAGS)  nbody.C 
nbody : nbody.C  nbody_particle.h $(LIBOBJS)
	g++  -DEVOLVE    $(CPPFLAGS) -o nbody nbody.C $(LIBOBJS) $(PGLIB)
nbody.dfs : nbody.C  nbody_particle.h BHtree.h nbody_system.h BHtree.C $(OTHERSRCS) $(OTHEROBJS) 
	g++  -DEVOLVE -DDFS_TREE_WALK   $(CPPFLAGS) -o nbody.dfs nbody.C  $(OTHERSRCS) $(OTHEROBJS)   $(PGLIB)
nbody.dfs.a64fx : nbody.C  nbody_particle.h BHtree.h nbody_system.h BHtree.C $(OTHERSRCS) $(OTHEROBJS) 
	FCC -Nclang  -DEVOLVE -DDFS_TREE_WALK   $(CLANGFLAGS) -o nbody.dfs.a64fx nbody.C  $(OTHERSRCS) $(OTHEROBJS)   $(PGLIB)
nbody.st : nbody.C  nbody_particle.h BHtree.h nbody_system.h BHtree.C $(OTHERSRCS) $(OTHEROBJS) 
	g++  -DEVOLVE -DSIMD_TREE_WALK   $(CPPFLAGS) -o nbody.st nbody.C  $(OTHERSRCS) $(OTHEROBJS)   $(PGLIB)
nbody.stdfs : nbody.C  nbody_particle.h BHtree.h nbody_system.h BHtree.C $(OTHERSRCS) $(OTHEROBJS) 
	g++  -DEVOLVE -DSIMD_DFS_TREE_WALK   $(CPPFLAGS) -o nbody.stdfs nbody.C  $(OTHERSRCS) $(OTHEROBJS)   $(PGLIB)
nbody.st.a64fx  : nbody.C  nbody_particle.h BHtree.h nbody_system.h BHtree.C $(OTHERSRCS) $(OTHEROBJS) 
	FCC -Nclang  -DEVOLVE -DSIMD_TREE_WALK   $(CLANGFLAGS) -o nbody.st.a64fx nbody.C  $(OTHERSRCS) $(OTHEROBJS)   $(PGLIB)
nbody.stdfs.a64fx  : nbody.C  nbody_particle.h BHtree.h nbody_system.h BHtree.C $(OTHERSRCS) $(OTHEROBJS) 
	FCC -Nclang  -DEVOLVE -DSIMD_DFS_TREE_WALK   $(CLANGFLAGS) -o nbody.stdfs.a64fx nbody.C  $(OTHERSRCS) $(OTHEROBJS)   $(PGLIB)
nbody.a64fx  : nbody.C  nbody_particle.h BHtree.h nbody_system.h BHtree.C $(OTHERSRCS) $(OTHEROBJS) 
	FCC -Nclang  -DEVOLVE   $(CLANGFLAGS) -o nbody.a64fx nbody.C  $(OTHERSRCS) $(OTHEROBJS)   $(PGLIB)
nbtest : nbtest.C  nbody_particle.h $(LIBOBJS)
	g++  -DEVOLVE    $(CPPFLAGS) -o nbtest nbtest.C $(LIBOBJS) $(PGLIB)
nbtest.clang : nbtest.C  nbody_particle.h $(OTHERSRCS)  $(OTHEROBJS)
	FCC -Nclang  -DEVOLVE    $(CLANGFLAGS) -o nbtest.clang nbtest.C  $(OTHERSRCS)  $(OTHEROBJS)
nbtest.a64fx  : nbtest.cc  nbody_particle.h $(OTHERSRCS)  $(OTHEROBJS) Makefile
	$(A64CCOMPILER)  -DEVOLVE    $(A64FLAGS) nbtest.cc   $(OTHERSRCS)  $(OTHEROBJS)   -o nbtest.a64fx 
nbody_ext : nbody.C external.C gravity.C  nbody_particle.h $(LIBOBJSE)
	g++   -DEVOLVE -DEXTERNAL_FIELD   $(CPPFLAGS) -o nbody_ext nbody.C gravity.C external.C  $(LIBOBJSE) $(PGLIB)
nbody_g4 : nbody.C  nbody_particle.h $(LIBOBJSG4)
	g++  -DEVOLVE    $(CPPFLAGS) -o nbody_g4 nbody.C $(LIBOBJSG4) $(PGLIB) $(G4LIB)
nbody_g6 : nbody.C  nbody_particle.h $(LIBOBJSG6)
	g++  -DEVOLVE    $(CPPFLAGS) -o nbody_g6 nbody.C $(LIBOBJSG6) $(PGLIB) $(G6LIB)
nbody_h3 : nbody.C  nbody_particle.h $(LIBOBJSG4)
	g++  -DEVOLVE    $(CPPFLAGS) -o nbody_h3 nbody.C $(LIBOBJSG4) $(H3LIB)

BHtree.o : BHtree.C  nbody_particle.h BHtree.h 
	g++    $(CPPFLAGS) -c BHtree.C 
BHtree_g4.o : BHtree.C  nbody_particle.h BHtree.h 
	g++ -DHARP3   $(CPPFLAGS) -c BHtree.C
	mv BHtree.o  BHtree_g4.o 
BHtree_g6.o : BHtree.C  nbody_particle.h BHtree.h 
	g++ -DGRAPE6   $(CPPFLAGS) -c BHtree.C
	mv BHtree.o  BHtree_g6.o 
pgetopt.o : pgetopt.C  
	g++  $(CPPFLAGS)  -c pgetopt.C
samplelog:  nbody
	nbody -i samplein > samplelog
samplelog_g4:  nbody_g4
	nbody_g4 -i samplein > samplelog_g4

document: .document
	ls -al .document
.document: nbody_guide.tex
	latex nbody_guide
	latex2web.csh nbody_guide $(WEBDIR) 
	touch .document
.web: export.tar.gz .document
	/bin/cp -p export.tar.gz  $(WEBDIR)
	touch .web
web: .web
	ls -al .web

