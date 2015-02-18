OBJS = main.o

CC = g++
DEBUG = -g
LFLAGS = -Wall $(DEBUG)
CFLAGS = -Wall -c $(DEBUG)

vpath %.cpp ../Source
vpath %.hpp ../Source

all : projetPages

projetPages.exe : $(OBJS)
	$(CC) $(LFAGS) $(OBJS) -o projetPages.exe

%.o : %.cpp
	$(CC) $(CFLAGS) $<
