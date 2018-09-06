# Makefile
# Author     : Yueze Tan
# Email      : yut75@psu.edu
# Written for: CSE 597-002(Fall 2018), HW00
# Last modified @ 2018/08/22

NAME = HelloWorld

CC = g++
CFLAGS = -O2
LFLAGS = 

SRCPATH = ./src
OBJPATH = ./obj
BINPATH = ./bin

SOURCES = $(wildcard ${SRCPATH}/*.cpp)
OBJECTS = $(patsubst %.cpp, ${OBJPATH}/%.o, $(notdir ${SOURCES}))
BINNAME = $(NAME)

program: $(OBJECTS) | $(BINPATH)
	$(CC) $(OBJECTS) -o $(BINPATH)/$(BINNAME) $(LFLAGS)
	@echo ""
	@echo "Make completed."
	@echo ""

$(OBJPATH)/%.o: $(SRCPATH)/%.cpp | $(OBJPATH)
	$(CC) -c $< -o $@ $(CFLAGS)

$(BINPATH):
	@mkdir $(BINPATH)

$(OBJPATH):
	@mkdir $(OBJPATH)

clean:
	@rm -rf $(OBJPATH)
	@echo "Clean completed."
	@echo ""

