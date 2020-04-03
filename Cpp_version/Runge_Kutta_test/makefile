EXEC   = main_out

OPTIMIZE =  -O2  

OBJS   = main.o runge_kutta_solver.o

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

