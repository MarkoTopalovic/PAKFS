LIBDIR = /opt/paklibsi
MUMPS = $(LIBDIR)/MUMPS
include $(MUMPS)/Makefile.inc
# iz Makefile.inc
#FC = mpif90
#CC = icc

FLAGS= -O3 -I$(MUMPS)/include -I $(LIBDIR)/MPICHinstall/src/include -static

SRC:=./PAKS
SRCPP:=./cplusplus
DEPS = $(SRCPP)/*.h
MAKEd = PAKS/x64linux/

TARGETLIB = $(MAKEd)/libpak
TARGET = linuxexe/pak.exe

 
SUFFIX = for
SUFFIX90 = f90
SUFFIXC = c

temp := $(foreach dir,$(SRC),$(wildcard $(dir)/*.$(SUFFIX)))
temp90 := $(foreach dir,$(SRC),$(wildcard $(dir)/*.$(SUFFIX90)))
tempc := $(foreach dir,$(SRCPP),$(wildcard $(dir)/*.$(SUFFIXC)))

srcs := $(notdir $(temp))
srcs90 := $(notdir $(temp90))
srcsc := $(notdir $(tempc))

mak  := $(srcs:.for=.obj)
mak90  := $(srcs90:.f90=.obj)
makc  := $(srcsc:.c=.obj)

makd := $(addprefix $(MAKEd),$(mak))
makd90 := $(addprefix $(MAKEd),$(mak90))
makdc := $(addprefix $(MAKEd),$(makc))
RM = rm -f


LIB = $(TARGETLIB).a $(MUMPS)/lib/libdmumps.a $(MUMPS)/lib/libmumps_common.a $(LIBS) $(LIBDIR)/ATLAS/build/lib/libatlas.a -L/usr/lib64/ -lpthread $(MUMPS)/PORD/lib/libpord.a $(LIBDIR)/METIS/lib/libmetis.a

all: $(makd90) $(makdc) $(makd) makefile
	$(FC) $(FLAGS) $(LIB) -o $(TARGET)

$(MAKEd)%.obj: $(SRC)/%.$(SUFFIX90) # makefile
	$(FC) -c $(FLAGS)  $< -o $@
	ar rs $(TARGETLIB).a $@

$(MAKEd)%.obj: $(SRCPP)/%.$(SUFFIXC) $(DEPS) # makefile
	$(CC) -c $(FLAGS)  $< -o $@
	ar rs $(TARGETLIB).a $@

$(MAKEd)%.obj: $(SRC)/%.$(SUFFIX) # makefile
	$(FC) -c $(FLAGS)  $< -o $@
	ar rs $(TARGETLIB).a $@
	
clean:
	rm ./$(MAKEd)*.obj
	rm ./$(MAKEd)*.a
	rm ./*.mod
	rm ./linuxexe/pak.exe

# DO NOT DELETE THIS LINE -- make depend depends on it.
