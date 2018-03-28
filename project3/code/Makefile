## Compiler, tools and options
CC      = cc
# compiler flags: to avoid warning
CCFLAGS = -O3

## Files
OBJECTS = summa.o
TARGET  = summa
SRC = summa.c

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
