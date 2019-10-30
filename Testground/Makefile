# How to use 
# To make tests: make tests 
# To make all experiment binaries: make binaries

CXX = /usr/local/Cellar/gcc/9.1.0/bin/g++-9
CFLAGS = -O3 -std=c++11 -fopenmp

SRCS = Debug.cpp SequenceMinHash.cpp SignedRandomProjections.cpp L2Hash.cpp MinHash.cpp driver.cpp io.cpp MurmurHash.cpp reservoir.cpp util.cpp RACE.cpp 
SRCS_DIR = src/

BUILD_DIR = build/
BIN_DIR = bin/
INC := -I include

# List of target executables
TARGETS = main.cpp RACEThresholds.cpp RSReservoir.cpp KNNBuffer.cpp
TARGETS_DIR = targets/

# Everything beyond this point is determined from previous declarations, don't modify
OBJECTS = $(addprefix $(BUILD_DIR), $(SRCS:.cpp=.o))
BINARIES = $(addprefix $(BIN_DIR), $(TARGETS:.cpp=))

$(BUILD_DIR)%.o: $(SRCS_DIR)%.cpp
	$(CXX) $(INC) -c $(CFLAGS) $< -o $@

binaries: $(BINARIES)
targets: $(BINARIES)
exp: $(BINARIES)

$(BINARIES): $(addprefix $(TARGETS_DIR), $(TARGETS)) $(OBJECTS)
	$(CXX) $(INC) $(CFLAGS) $(OBJECTS) $(addsuffix .cpp,$(@:$(BIN_DIR)%=$(TARGETS_DIR)%)) -o $@

clean:
	rm -f $(OBJECTS); 
	rm -f $(BINARIES); 

.PHONY: clean targets binaries exp

