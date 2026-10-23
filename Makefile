# copyright: Michael Safyan
# modified by Falk Hildebrand

program_NAME := LCA
PROVE ?= prove
program_C_SRCS := $(wildcard src/*.c)
program_CXX_SRCS := $(wildcard src/*.cpp)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program_DEPS := $(program_OBJS:.o=.d)
program_INCLUDE_DIRS :=
program_LIBRARY_DIRS :=
program_LIBRARIES :=

CXXFLAGS += -Wall -Wextra -Wpedantic -O3 -std=c++20 -pthread
DEPFLAGS := -MMD -MP
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDLIBS += $(foreach library,$(program_LIBRARIES),-l$(library)) -lz

.PHONY: all check clean distclean

all: $(program_NAME)

check: $(program_NAME) tests/regression
	./tests/regression
	LCA_BIN=./$(program_NAME) $(PROVE) tests/cli_regression.t

# The unit tests link every source but the program's main().
test_SRCS := $(filter-out src/LCA.cpp,$(program_CXX_SRCS))

tests/regression: tests/regression.cpp $(test_SRCS) $(wildcard src/*.h)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) tests/regression.cpp $(test_SRCS) $(LDFLAGS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(program_NAME): $(program_OBJS)
	$(LINK.cc) $(program_OBJS) $(LDLIBS) -o $(program_NAME)

-include $(program_DEPS)

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) tests/regression tests/regression.exe
	@- $(RM) $(program_OBJS)
	@- $(RM) $(program_DEPS)

distclean: clean
