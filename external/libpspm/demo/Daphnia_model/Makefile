#-------------------------------------------------------------------------------
# executable name
FILE = 
TARGET = $(patsubst %.cpp, %.exec, $(FILE)) 

# files
# SRCFILES  :=  $(wildcard ../../src/*.cpp plant/*.cpp)
SRCFILES := $(FILE)
HEADERS := #$(wildcard ../../src/*.tpp) $(wildcard ../../include/*.h) 
# ------------------------------------------------------------------------------

# paths
#CUDA_INSTALL_PATH ?= /usr/local/cuda#-5.0

# include and lib dirs (esp for cuda)
INC_PATH := -I../../include 
LIB_PATH := -L../../lib 

# flags
CPPFLAGS = -O3 -g -pg -std=c++17 -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable
LDFLAGS =  -g -pg

# libs
LIBS = -lpspm	 #-lgsl -lgslcblas 	# additional libs
#LIBS = -lcudart 			# cuda libs 		

# files
OBJECTS = $(patsubst %.cpp, %.o, $(SRCFILES))


all: dir $(TARGET)	

dir: 
	mkdir -p build 

$(TARGET): $(OBJECTS)	
	g++ $(LDFLAGS) -o $(TARGET) $(LIB_PATH) $(OBJECTS) $(LIBS) 

$(OBJECTS): %.o : %.cpp $(HEADERS)
	g++ -c $(CPPFLAGS) $(INC_PATH) $< -o $@ 

clean:
	rm -f *.o *.exec gmon.out 

dataclean:	
	rm -f *.txt
	
re: clean all









#-gencode=arch=compute_10,code=\"sm_10,compute_10\"  -gencode=arch=compute_20,code=\"sm_20,compute_20\"  -gencode=arch=compute_30,code=\"sm_30,compute_30\" 

#-W -Wall -Wimplicit -Wswitch -Wformat -Wchar-subscripts -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wno-unused-function 
#-m64 -fno-strict-aliasing 
#-I. -I/usr/local/cuda/include -I../../common/inc -I../../../shared//inc 
#-DUNIX -O2


#g++ -fPIC   -m64 -o ../../bin/linux/release/swarming_chasing_predator obj/x86_64/release/genmtrand.cpp.o  obj/x86_64/release/simpleGL.cu.o  -L/usr/local/cuda/lib64 -L../../lib -L../../common/lib/linux -L../../../shared//lib -lcudart   
#-lGL -lGLU -lX11 -lXi -lXmu -lGLEW_x86_64 -L/usr/X11R6/lib64 -lGLEW_x86_64 -L/usr/X11R6/lib64 -lglut 
#-L/usr/local/cuda/lib64 -L../../lib -L../../common/lib/linux -L../../../shared//lib -lcudart 
#-L/usr/lib -lgsl -lgslcblas 
#-lcutil_x86_64  -lshrutil_x86_64 




#CXXWARN_FLAGS := \
#	-W -Wall \
#	-Wimplicit \
#	-Wswitch \
#	-Wformat \
#	-Wchar-subscripts \
#	-Wparentheses \
#	-Wmultichar \
#	-Wtrigraphs \
#	-Wpointer-arith \
#	-Wcast-align \
#	-Wreturn-type \
#	-Wno-unused-function \
#	$(SPACE)

#CWARN_FLAGS := $(CXXWARN_FLAGS) \
#	-Wstrict-prototypes \
#	-Wmissing-prototypes \
#	-Wmissing-declarations \
#	-Wnested-externs \
#	-Wmain \
#	
	
#HEADERS  := $(wildcard *.h)
	
