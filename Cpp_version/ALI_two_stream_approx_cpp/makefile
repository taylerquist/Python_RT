EXEC   = main_out

OPTIMIZE =  -O2  

OBJS   = main.o tau_fn.o quad_int.o matrix_build.o s_interp.o tri_solver.o mat_mult.o

.SUFFIXES : .c .cpp .o

INCL   = -I./
LIBS   = -lm 


FLAGS = 
CFLAGS 	  = $(OPTIMIZE) $(FLAGS) -m64
CXXFLAGS  = $(OPTIMIZE) $(FLAGS) -m64
LDFLAGS	  = -m64


%.o:	%.c
		$(CC) $(CFLAGS)  $(INCL)  -c $< -o $@ 

%.o:	%.cpp
		$(CXX) $(CXXFLAGS)  $(INCL)  -c $< -o $@ 

$(EXEC): $(OBJS) 
	 	 $(CXX) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXEC) $(INCL)   

#$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

