# Compiler and flags
CXX = gcc
CXXFLAGS = -Wall

# Object files
OBJECTS = sequenceAVLTree.o filterSequence.o

# Target executable
sequenceAVLTree: $(OBJECTS)
	$(CXX) -g $(CXXFLAGS) -o sequenceAVLTree $(OBJECTS)

# Object file generation for sequenceAVLTree.c
sequenceAVLTree.o: sequenceAVLTree.c sequenceAVLTree.h
	$(CXX) -g $(CXXFLAGS) -c sequenceAVLTree.c

# Object file generation for filterSequence.c
filterSequence.o: filterSequence.c sequenceAVLTree.h
	$(CXX) -g $(CXXFLAGS) -c filterSequence.c

# Clean command to remove object files and the executable
clean:
	rm -f *.o
	rm -f sequenceAVLTree
