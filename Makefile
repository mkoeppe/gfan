ALL: default

# When -g is used, it seems that stackdump does not work.
# GPROFFLAG = -pg
GPROFFLAG = -g
#GPROFFLAG = 

PREFIX?=/usr/local

BINDIR=$(PREFIX)/bin

RELEASEDIR = gfan0.7
MAIN       = gfan

GCATSPATH   = ./


ifeq ($(stackdump),)
STACTDUMP_OPTIONS =
else
STACTDUMP_OPTIONS = -DSTACKDUMP_ENABLED -D__assert_fail=__assert_fail2
endif


ifeq ($(sagepath),)
SAGE_LINKOPTIONS = 
SAGE_INCLUDEOPTIONS =
else
SAGE_LINKOPTIONS = -L$(sagepath)/ -lpython2.6 -lcsage -lsingular
SAGE_INCLUDEOPTIONS = -I $(sagepath)/
SAGE_OBJECTS = sage.o sage_link.so

sage_link.so: sage_link.pyx setup.py
	python setup.py build_ext --inplace --pyrex-include-dirs=$(SAGE_ROOT)/devel/sage/
endif

ifeq ($(gmppath),)
GMP_LINKOPTIONS = -lgmp
GMP_INCLUDEOPTIONS =
else
GMP_LINKOPTIONS = $(gmppath)/lib/libgmp.a
GMP_INCLUDEOPTIONS = -I $(gmppath)/include
endif

ifeq ($(cddnoprefix),)
CDDDEFINE_PREFIX =
else
CDDDEFINE_PREFIX = -DNOCDDPREFIX
endif

ifeq ($(cddpath),)
CDD_LINKOPTIONS = -L/usr/local -lcddgmp
CDD_INCLUDEOPTIONS =
else
CDD_LINKOPTIONS = $(cddpath)/lib/libcddgmp.a
CDD_INCLUDEOPTIONS = -I $(cddpath)/include
endif

ifeq ($(soplex),)
SOPLEX_PATH =
SOPLEX_LINKOPTIONS =
SOPLEX_INCLUDEOPTIONS =
SOPLEX_OBJECTS =
else
SOPLEX_PATH = $(HOME)/math/software/soplex-1.3.2
SOPLEX_LINKOPTIONS = $(SOPLEX_PATH)/lib/libsoplex.linux.x86_64.gnu.opt.a -lz
#SOPLEX_LINKOPTIONS = -lz $(SOPLEX_PATH)/lib/libsoplex.linux.x86.gnu.opt.a
#SOPLEX_LINKOPTIONS = -lz $(SOPLEX_PATH)/lib/libsoplex.darwin.x86.gnu.opt.a
SOPLEX_INCLUDEOPTIONS = -I $(SOPLEX_PATH)/src
SOPLEX_OBJECTS = lp_soplexcdd.o
endif

ifeq ($(vectorize),)
VECTORISE =
else
VECTORISE = -mavx -msse2
endif

# rememember to adjust USEFACTORY in field_rationalfunctions2
ifeq ($(singular),)
ifeq ($(factory),)
SINGULAR_PATH =
SINGULAR_LINKOPTIONS =
SINGULAR_INCLUDEOPTIONS =
SINGULAR_OBJECTS = src/polynomialgcd.o 
else
SINGULAR_PATH =
SINGULAR_LINKOPTIONS = -lcf -lcfmem
SINGULAR_INCLUDEOPTIONS = -I /usr/local/include/factory/
SINGULAR_OBJECTS = src/polynomialgcd.o src/ftmpl_inst.o
endif
else
#SINGULAR_PATH = $(HOME)/math/software/Singular-3-1-0
#SINGULAR_LINKOPTIONS =  -L$(SINGULAR_PATH)/Singular -lsingular -lncurses -lreadline
#SINGULAR_INCLUDEOPTIONS = -I $(SINGULAR_PATH)/kernel -I $(SINGULAR_PATH)/omalloc
#SINGULAR_OBJECTS = src/singular.o src/singularconversion.o
SINGULAR_PATH = $(HOME)/math/software/Singular-svn/trunk/x86_64-Linux
SINGULAR_LINKOPTIONS =  -L$(SINGULAR_PATH)/lib -lsingular -lncurses -lreadline  -lcf -lcfmem
SINGULAR_INCLUDEOPTIONS = -I $(SINGULAR_PATH)/include -I $(SINGULAR_PATH)/include/omalloc
SINGULAR_OBJECTS = src/ftmpl_inst.o src/singular.o src/singularconversion.o
#Run the following line before running gfan
#export LD_LIBRARY_PATH="/home/anders/math/software/Singular-svn/trunk/x86_64-Linux/lib/:${LD_LIBRARY_PATH}"
endif

# To produce factory templates:
#g++ -c /home/anders/math/software/factory-3-1-7/ftmpl_inst.cc  -fno-implicit-templates -I /usr/local/include/factory/ -I/home/anders/math/software/factory-3-1-7/ -O2 -fomit-frame-pointer -o ftmpl_inst.o

ADDITIONALLINKOPTIONS = $(CDD_LINKOPTIONS) $(GMP_LINKOPTIONS) $(SOPLEX_LINKOPTIONS) $(SINGULAR_LINKOPTIONS) $(SAGE_LINKOPTIONS) 
ADDITIONALINCLUDEOPTIONS = $(CDD_INCLUDEOPTIONS) $(GMP_INCLUDEOPTIONS) $(SOPLEX_INCLUDEOPTIONS) $(SINGULAR_INCLUDEOPTIONS) $(SAGE_INCLUDEOPTIONS)


MKDIR=mkdir -p

# PREFIX = /usr/local/gcc/6.2.0/bin/
# PREFIX = /usr/local/gcc-8.1/bin/
PREFIX =
SHELL       = /bin/sh
#ARCH        = LINUX

CC          = $(PREFIX)gcc
CLINKER     = $(CC)
CXX         = $(PREFIX)g++
CCLINKER    = $(CXX)

#CC          = $(PREFIX)gcc-8.1
#CLINKER     = $(CC)
#CXX         = $(PREFIX)g++-8.1
#CCLINKER    = $(CXX)

#OPTFLAGS    = -O2 -DGMPRATIONAL -DNDEBUG
# Note that gcc produces wrong code with -O3
# OPTFLAGS    =  -DGMPRATIONAL -Wuninitialized -O3 -fno-omit-frame-pointer -fdelete-null-pointer-checks -fno-inline-functions -fno-inline-small-functions#-fno-common  #-fno-inline-functions #-D_GLIBCXX_DEBUG #-O2	 #-O3 -fno-guess-branch-probability #-DNDEBUG
OPTFLAGS    =  -DGMPRATIONAL -Wuninitialized -O3  $(VECTORISE)  -finline-limit=1000 -ffast-math -fno-omit-frame-pointer -fno-inline-functions -fdelete-null-pointer-checks #-D_GLIBCXX_DEBUG #-fno-inline-functions -fno-inline-small-functions#-fno-common  #-fno-inline-functions #-D_GLIBCXX_DEBUG #-O2	 #-O3 -fno-guess-branch-probability #-DNDEBUG
#OPTFLAGS    =  -DGMPRATIONAL -Wuninitialized -fno-omit-frame-pointer -O3 -mavx -msse2  -finline-limit=1000 -ffast-math -Wuninitialized # -fno-guess-branch-probability #-DNDEBUG -ftree-vectorizer-verbose=2
#OPTFLAGS    =  -DGMPRATIONAL -Wuninitialized -fno-omit-frame-pointer -O1             -fno-guess-branch-probability
 #-DNDEBUG
#OPTFLAGS    =  -DGMPRATIONAL -Wuninitialized -fno-omit-frame-pointer -O3 -mavx -msse2 -ftree-vectorizer-verbose=2 -finline-limit=1000 -ffast-math #-DNDEBUG
#OPTFLAGS    =  -DGMPRATIONAL -Wuninitialized -fno-omit-frame-pointer -O3 -mavx -msse2 -ftree-vectorizer-verbose=2 -march=native -unroll-loops --param max-unroll-times=4 -ffast-math #-DNDEBUG
#OPTFLAGS    =  -DGMPRATIONAL -Wuninitialized -fno-omit-frame-pointer -O3 -msse2 -ftree-vectorizer-verbose=2 -ffast-math #-DNDEBUG
#OPTFLAGS    =  -DGMPRATIONAL -Wuninitialized -fno-omit-frame-pointer -O3 -mavx -msse2 -ftree-vectorizer-verbose=2 -ffast-math #-DNDEBUG

CFLAGS	  = $(OPTFLAGS) $(GPROFFLAG) $(STACTDUMP_OPTIONS) $(ADDITIONALINCLUDEOPTIONS) -g $(CDDDEFINE_PREFIX) -std=c++2a #-pedantic
#CFLAGS	  = $(OPTFLAGS) $(GPROFFLAG) $(STACTDUMP_OPTIONS) $(ADDITIONALINCLUDEOPTIONS) -std=c++0x -g $(CDDDEFINE_PREFIX) #-pedantic
#CFLAGS	  = $(OPTFLAGS) $(GPROFFLAG) $(STACTDUMP_OPTIONS) $(ADDITIONALINCLUDEOPTIONS) -D_GLIBCXX_DEBUG -g $(CDDDEFINE_PREFIX) #-pedantic
CCFLAGS	  = $(CFLAGS)
FFLAGS	  = $(OPTFLAGS)

CATSOBJECTS =	lp_cdd.o \
		parser.o \
		field.o \
		monomial.o \
		printer.o	\
		polynomial.o \
		termorder.o \
		term.o \
		vektor.o \
		division.o \
		buchberger.o \
		wallideal.o \
		lp.o \
		enumeration.o \
		ep_standard.o \
		ep_xfig.o \
		reversesearch.o \
		application.o \
		timer.o \
		renderer.o \
		field_rationals.o \
		symmetry.o \
		breadthfirstsearch.o \
		genericwalk.o \
		minkowskisum.o \
		newtonpolytope.o \
		tropical.o \
		dimension.o \
		bergman.o \
		subspace.o \
		polyhedralcone.o \
		gfanapplication.o \
		polyhedralfan.o \
		tropical2.o \
		field_zmodpz.o \
		binomial.o \
		matrix.o \
		latticeideal.o \
		scarf.o \
		xfig.o \
		halfopencone.o \
		lll.o \
		multiplicity.o \
		substitute.o \
		polymakefile.o \
		saturation.o \
		determinant.o \
		polynomialring.o \
		log.o \
		tropicalbasis.o \
		symmetriccomplex.o \
		linalg.o \
		minors.o \
		continuedfractions.o \
		triangulation.o \
		minkowskidual.o \
		regularsubdivision.o \
		fieldlp.o \
		field_rationalfunctions.o \
		tropical_weildivisor.o \
		intsinpolytope.o\
		lattice.o \
		graph.o \
		restrictedautoreduction.o \
		tropicaltraverse.o \
		groebnerengine.o \
		ge_gfan.o \
		nbody.o \
		codimoneconnectedness.o \
		tropicalmap.o \
		traverser_tropical.o \
		traverser_groebnerfan.o \
		field_rationalfunctions2.o \
	mixedvolume.o \
	traverser_stableintersection.o \
	traverser_secondaryfan.o \
	linalgfloat.o \
	primarydecomposition.o \
	tropicaldeterminant.o \
	determinantpoly.o \
	traverser_sphere.o \
	gfanlib_zcone.o \
	gfanlib_symmetry.o \
	gfanlib_symmetriccomplex.o \
	gfanlib_polyhedralfan.o \
	gfanlib_zfan.o \
	gfanlib_polymakefile.o \
	gfanlib_mixedvolume.o \
	gfanlib_circuittableint.o \
	gfanlib_paralleltraverser.o \
	padic.o \
	integergb.o \
	traverser_resultantfan.o \
	bsptree.o \
	traverser_resultantfanspecialization.o \
	myassert.o \
	traverser_bsptree.o \
	gfanlib_traversal.o \
	tropicalcurve.o \
	packedmonomial.o \
	gmpallocator.o \
	gfanlib_memoryresource.o \
	gfanlib_hypersurfaceintersection.o \
	divisionobject.o \
	gfanlibglue.o \
#	gfanlib_circuittableinteger.o \
#		symmetrictraversal.o \
	#	restrictedgfan.o \

APPDELETEOBJECTS = 		app_add.o \
		app_berndssuggestion.o \
		app_grassmanndata2.o \
		app_grassmanndata3.o \
		app_construction.o \
		app_checkridges.o \
		app_edwinsconjecture.o \
		app_fvector.o \
		app_grassmanndata.o \
		app_groupfacetbinomials.o \
		app_istriangulation.o \
		app_latticetest.o \
		app_markpolynomialset.o \
		app_moeckel.o \
		app_polytopetopolynomial.o \
		app_rendernewtonpolytope.o \
		app_tropical.o \
		app_xfigconstruction.o \
		app_liststandardmonomials.o \
#		app_isrefinement.o \    # needs to be fixed so that it compiles with gcc version 2.96 (legolas.imf.au.dk)


APPOBJECTS = app_main.o \
		app_buchberger.o \
		app_doesidealcontain.o \
		app_facets.o \
		app_groebnercone.o \
		app_homogeneityspace.o \
		app_homogenize.o \
		app_initialforms.o \
		app_interactive.o \
		app_isgroebnerbasis.o \
		app_ismarkedgroebnerbasis.o \
		app_krulldimension.o \
		app_leadingterms.o \
		app_multiplymatrix.o \
		app_polynomialsetunion.o \
		app_render.o \
		app_renderstaircase.o \
		app_stats.o \
		app_substitute.o \
		app_supportindices.o \
		app_tolatex.o \
		app_transposematrix.o \
		app_tropicalbasis.o \
		app_tropicalintersection.o \
		app_tropicalstartingcone.o \
		app_tropicaltraverse.o \
		app_walk.o \
		app_weightvector.o \
		app_scarfisgeneric.o \
		app_scarfvisualize.o \
		app_scarfcomplex.o \
		app_sturmsequence.o \
		app_latticeideal.o \
		app_lll.o \
		app_tropicalmultiplicity.o \
		app_idealintersection.o \
		app_test.o \
		app_saturation.o \
		app_idealproduct.o \
		app_representatives.o \
		app_tropicallifting.o \
		app_topolyhedralfan.o \
		app_tropicalbruteforce.o \
		app_secondaryfan.o \
		app_composepermutations.o \
		app_minors.o \
		app_tropicalrank.o \
		app_minkowski.o \
		app_triangulate.o \
		app_tropicallinearspace.o \
		app_combinerays.o \
		app_regularsubdivision.o \
		app_lpsolve.o \
		app_tropicalweildivisor.o \
		app_lattice.o \
		app_intsinpolytope.o\
		app_tropicalevaluation.o \
		app_smalessixth.o \
		app_smalessixth2.o \
		app_nbody.o \
		app_spolynomial.o \
		app_link.o \
		app_normalfancleanup.o \
		app_tropicalfunction.o \
		app_volume.o \
		app_isconnected.o \
		app_tropicalhypersurface.o \
		app_product.o \
		app_commonrefinement.o \
		app_tropicalimage.o \
		app_groebnerfan.o \
		app_fanhomology.o \
		app_genericlinearchange.o \
		app_mixedvolume.o \
		app_fiberpolytope.o \
		app_symmetries.o \
		app_evaluate.o \
		app_exponentlattice.o \
		app_minimalassociatedprimes.o \
		app_realroots.o \
		app_initialdeterminant.o \
		app_fansubfan.o \
		app_fancones.o \
		app_issmooth.o \
		app_fancoarsening.o \
		app_pointconfiguration.o \
		app_librarytest.o \
		app_padic.o \
		app_integergb.o \
		app_matrixproduct.o \
		app_traversetropicalintersection.o \
		app_markpolynomialset.o \
		app_tropicalhypersurfacereconstruction.o \
		app_resultantfan.o \
		app_isbalanced.o \
		app_polytopealgebra.o \
		app_debug.o \
		app_randompolynomials.o \
		app_tropicalcurve.o \
		app_tropicalhomotopy.o \
		app_integerfactorization.o \
		app_tropicalvarietyspan.o \
		app_chowbetti.o \
		app_anton.o \
		app_components.o \
		app_tropicalprevarietycomponents.o \
		app_tropicalprevariety.o \
		app_anders.o \

GFANLIBFILES= 	gfanlib.h \
				gfanlib_polyhedralfan.cpp \
				gfanlib_q.h \
				gfanlib_symmetry.h \
				gfanlib_zcone.cpp \
				gfanlib_z.h \
				gfanlib_matrix.h \
				gfanlib_polyhedralfan.h \
				gfanlib_symmetriccomplex.cpp \
				gfanlib_zcone.h \
				gfanlib_polymakefile.cpp \
				gfanlib_symmetriccomplex.h \
				gfanlib_zfan.cpp \
				gfanlib_polymakefile.h \
				gfanlib_symmetry.cpp \
				gfanlib_vector.h \
				gfanlib_zfan.h \
				gfanlib_circuittableint.h \
				gfanlib_circuittableint.cpp \
				gfanlib_mixedvolume.h \
				gfanlib_mixedvolume.cpp \
				gfanlib_paralleltraverser.h \
				gfanlib_paralleltraverser.cpp \
				gfanlib_tropicalhomotopy.h \
				gfanlib_memoryresource.h \
				gfanlib_memoryresource.cpp \
	#			gfanlib_tableau.h \   we do not yet support this one


EXECS	  = $(MAIN)

# When compiling symmetrictraversal.cpp as any other file, it causes the program to crash in certain situations.
# (compiling with gcc version 4.7.2 and running gfan _tropicaltraverse on a starting cone for Grassmann3_7)
# Either this is a bug in the code or in the compiler. The bug disappears by compiling with -fno-guess-branch-probability
src/symmetrictraversal.o: src/symmetrictraversal.cpp
	$(CXX) $(CFLAGS) -fno-guess-branch-probability  -c src/symmetrictraversal.cpp -o src/symmetrictraversal.o
# If compiling with clang, use the line below instead:
#	$(CXX) $(CFLAGS) -c src/symmetrictraversal.cpp -o src/symmetrictraversal.o

# Define suffixes to make the program compile on legolas.imf.au.dk :
.SUFFIXES: .o .cpp .c

OBJECTS1 = 	$(addprefix src/,$(SOPLEX_OBJECTS)) \
		$(SINGULAR_OBJECTS) \
		$(SAGE_OBJECTS) \
		$(addprefix src/,$(CATSOBJECTS)) \
		$(addprefix src/,$(APPOBJECTS)) \
		src/symmetrictraversal.o

OBJECTS2 = src/app_test.o \
		   src/gfanapplication.o \
		   src/application.o \
		   src/polynomial.o \
		   src/polynomialring.o \
		   src/term.o \
		   src/termorder.o \
		   src/monomial.o \
		   src/printer.o \
		   src/field.o \
		   src/vektor.o \
		   src/log.o \
		   src/gfanlib_polymakefile.o \
		   src/linalg.o \
		   src/field_rationals.o \
		   src/linalgfloat.o \
		   src/myassert.o \
		   src/parser.o \
		   src/lp.o \
		   src/timer.o \
		   src/lp_cdd.o \
		   src/gfanlib_circuittableint.o \
		   src/gfanlib_paralleltraverser.o \
		   src/gfanlib_memoryresource.o \
		   src/symmetry.o \
		   src/division.o \
		   src/polyhedralcone.o \
		   src/subspace.o \
		   src/matrix.o \
		   src/polymakefile.o \
		   src/polyhedralfan.o \
		   src/symmetriccomplex.o \
		   src/symmetrictraversal.o \
		   src/wallideal.o \
		   src/enumeration.o \
		   src/buchberger.o \
		   src/traverser_groebnerfan.o \
		src/ge_gfan.o \
		src/tropical.o \
		src/tropical2.o \
		src/field_zmodpz.o \
		src/continuedfractions.o \
		src/determinant.o \
		src/reversesearch.o \
		src/groebnerengine.o \
		src/graph.o \
		src/dimension.o \
		src/polynomialgcd.o \
		src/field_rationalfunctions2.o \
		src/halfopencone.o \
		src/binomial.o \
		src/newtonpolytope.o \
		src/saturation.o \
		src/tropicalcurve.o \
		src/triangulation.o \
		src/breadthfirstsearch.o \
		src/gfanlib_hypersurfaceintersection.o \
		src/codimoneconnectedness.o


OBJECTS = $(OBJECTS1)

all: $(MAIN)

$(BINDIR): $(PREFIX)
	$(MKDIR) $(BINDIR)

$(PREFIX):
	$(MKDIR) $(PREFIX)

default: $(OBJECTS) $(ADDITIONALOBJECTS) $(EXECS)

$(MAIN): $(OBJECTS)
#	$(CCLINKER) $(OBJECTS) $(ADDITIONALLINKOPTIONS) $(GPROFFLAG) -lpthread  -o $(MAIN)
	$(CCLINKER) $(OBJECTS) $(ADDITIONALLINKOPTIONS) $(GPROFFLAG) -lpthread -rdynamic -o $(MAIN)

release:
	rm -f -r $(RELEASEDIR)/*
	mkdir -p $(RELEASEDIR)
	mkdir -p $(RELEASEDIR)/src
	cp src/*.cpp $(RELEASEDIR)/src
	cp src/*.h $(RELEASEDIR)/src
	cp Makefile $(RELEASEDIR)
#	cp macmake $(RELEASEDIR)
	cp README $(RELEASEDIR)
	cp LICENSE $(RELEASEDIR)
	cp COPYING $(RELEASEDIR)

	mkdir -p $(RELEASEDIR)/gfanlib
#	cp gfanlib/Makefile $(RELEASEDIR)/gfanlib/
	cp gfanlib/Makefile.in $(RELEASEDIR)/gfanlib/
	cp gfanlib/README.txt $(RELEASEDIR)/gfanlib/
	cp gfanlib/configure $(RELEASEDIR)/gfanlib/
	cp gfanlib/configure.ac $(RELEASEDIR)/gfanlib/

	cp Doxyfile $(RELEASEDIR)
	echo const char *GFAN_RELEASEDIR=\"$(RELEASEDIR)\"";" const char *GFAN_FORKTIME= >$(RELEASEDIR)/src/versioninfo.h
	date "+\"%s %a %h %d %H:%M:%S %Y\"" >>$(RELEASEDIR)/src/versioninfo.h
	echo ";" >>$(RELEASEDIR)/src/versioninfo.h
	mkdir -p $(RELEASEDIR)/examples/
#General examples:
	cp examples/2x2of2x3 $(RELEASEDIR)/examples/
	cp examples/2x2of2x4 $(RELEASEDIR)/examples/
	cp examples/2x2of3x3 $(RELEASEDIR)/examples/
	cp examples/2x2of4x4 $(RELEASEDIR)/examples/
	cp examples/4x4of4x5 $(RELEASEDIR)/examples/
	cp examples/4x4of5x5 $(RELEASEDIR)/examples/
	cp examples/6x6-subPfaffians $(RELEASEDIR)/examples/
	cp examples/cyclic4 $(RELEASEDIR)/examples/
	cp examples/linhyper5_2 $(RELEASEDIR)/examples/
	cp examples/linhyper5_2.cone $(RELEASEDIR)/examples/
	cp examples/pablo $(RELEASEDIR)/examples/
	cp examples/symmetryTest $(RELEASEDIR)/examples/
#Examples in Groebner fan paper:
	cp examples/examplePaper $(RELEASEDIR)/examples/
	cp examples/sturmfels3.9 $(RELEASEDIR)/examples/
	cp examples/3x3of3x4 $(RELEASEDIR)/examples/
	cp examples/3x3of3x5 $(RELEASEDIR)/examples/
	cp examples/3x3of4x4 $(RELEASEDIR)/examples/
#	cp examples/3x3of4x4sym $(RELEASEDIR)/examples/
	cp examples/grassmann2_5 $(RELEASEDIR)/examples/
	cp examples/cyclic5 $(RELEASEDIR)/examples/
#	cp examples/J4 $(RELEASEDIR)/examples/
#Examples useful for tropical computations:
#	cp examples/grassmann2_5 $(RELEASEDIR)/examples/
	cp examples/grassmann2_5.cone $(RELEASEDIR)/examples/
	cp examples/grassmann2_6 $(RELEASEDIR)/examples/
	cp examples/grassmann2_6.cone $(RELEASEDIR)/examples/
	cp examples/grassmann3_6 $(RELEASEDIR)/examples/
	cp examples/grassmann3_6.cone $(RELEASEDIR)/examples/
#Examples in tropical paper:
	cp examples/hankel3x3of4x4 $(RELEASEDIR)/examples/
	cp examples/hankel3x3of4x4.cone $(RELEASEDIR)/examples/
	cp examples/hankel3x3of4x5 $(RELEASEDIR)/examples/
	cp examples/hankel3x3of4x5.cone $(RELEASEDIR)/examples/
#	cp examples/hankel3x3of5x5 $(RELEASEDIR)/examples/
#	cp examples/hankel3x3of5x5.cone $(RELEASEDIR)/examples/
	cp examples/3x3of3x5.cone $(RELEASEDIR)/examples/
#	cp examples/3x3of4x4sym $(RELEASEDIR)/examples/
	cp examples/3x3of4x4sym.cone $(RELEASEDIR)/examples/
#	cp examples/3x3of5x5sym $(RELEASEDIR)/examples/
#	cp examples/3x3of5x5sym.cone $(RELEASEDIR)/examples/
	cp examples/commat2x2 $(RELEASEDIR)/examples/
	cp examples/commat2x2.cone $(RELEASEDIR)/examples/
#	cp examples/commat3x3 $(RELEASEDIR)/examples/
#	cp examples/commat3x3.cone $(RELEASEDIR)/examples/

#	cp -r testsuite $(RELEASEDIR)/          # The line below is equivalent but does not copy hidden files.
	rsync -av --exclude=".*" testsuite $(RELEASEDIR)/

	mkdir -p $(RELEASEDIR)/doc/
	cp doc/Makefile $(RELEASEDIR)/doc/
	cp doc/*.bib $(RELEASEDIR)/doc/
	cp doc/*.bbl $(RELEASEDIR)/doc/
	cp doc/*.tex $(RELEASEDIR)/doc/
	cp doc/*.dvi $(RELEASEDIR)/doc/
	cp doc/*.eps $(RELEASEDIR)/doc/
	cp doc/*.bst $(RELEASEDIR)/doc/
	cp doc/Makefile $(RELEASEDIR)/doc/
	mkdir -p $(RELEASEDIR)/homepage/
	cp webpage/*.png $(RELEASEDIR)/homepage/
	cp webpage/*.html $(RELEASEDIR)/homepage/
	cp webpage/Makefile $(RELEASEDIR)/homepage/
	mkdir -p $(RELEASEDIR)/homepage/presentation
	cp webpage/presentation/*.fig $(RELEASEDIR)/homepage/presentation/
	cp webpage/presentation/*.tex $(RELEASEDIR)/homepage/presentation/
	cp webpage/presentation/*.eps $(RELEASEDIR)/homepage/presentation/
	cp webpage/presentation/*.bib $(RELEASEDIR)/homepage/presentation/
	cp webpage/presentation/*.bbl $(RELEASEDIR)/homepage/presentation/
	cp webpage/presentation/*.ps $(RELEASEDIR)/homepage/presentation/
	cp webpage/presentation/*.pdf $(RELEASEDIR)/homepage/presentation/
	cp webpage/presentation/Makefile $(RELEASEDIR)/homepage/presentation/

	tar -c $(RELEASEDIR) > $(RELEASEDIR).tar  
	gzip $(RELEASEDIR).tar

clean:
	/bin/rm -f src/*.o $(EXECS) $(MAIN)
install: $(BINDIR)
	cp $(EXECS) $(BINDIR)
#	cp $(EXECS) /usr/local/bin
#	./gfan installlinks --path $(BINDIR)/
	cd $(BINDIR) && ./gfan installlinks
library:
	cd gfanlib && \
	/bin/rm -f *.a && \
	/bin/rm -f *.o && \
	/bin/rm -f *~ && \
	/bin/rm -f a.out && \
	/bin/rm -f Makefile && \
	/bin/rm -f config.status && \
	/bin/rm -f config.log
	cd ..
	cp $(addprefix src/,$(GFANLIBFILES)) gfanlib/
	tar zcf -  gfanlib --exclude=".*" > gfanlib0.6.2.tar.gz
check:
	./gfan _test
.c.o:
	$(CC) $(CFLAGS) -c $< -o $(patsubst %.c,%.o,$<)
.cpp.o:
	$(CXX) $(CFLAGS) -c $< -o $(patsubst %.cpp,%.o,$<)

# wget http://ftp.sunet.se/pub/gnu/gmp/gmp-4.2.2.tar.gz
# tar -xzvf gmp-4.2.2.tar.gz
# cd gmp-4.2.2
# ./configure --prefix=$HOME/gmp
# make
# make install
# make check
# cd..

# wget ftp://ftp.ifor.math.ethz.ch/pub/fukuda/cdd/cddlib-094d.tar.gz  # THIS LINE DOES NOT WORK!
# tar -xzvf cddlib-094d.tar.gz
# cd cddlib-094d
# ./configure --prefix="$HOME/cddlib" CFLAGS="-I$HOME/gmp/include -L$HOME/gmp/lib"
