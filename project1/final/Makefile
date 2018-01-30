## Compiler, tools and options
CC      = cc
# compiler flags: to avoid warning
CCFLAGS = -hmsglevel_4

## Files
OBJECTS = matmul.o
TARGET  = matmul
SRC = matmul.c papi_timer.c

## Implicit rules
.SUFFIXES: .c

#.c.o:
#	$(CC) -c $(CCFLAGS) $<

## Build rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CCFLAGS) -o $(OBJECTS) $(SRC)

clean:
	rm -f $(OBJECTS) 
