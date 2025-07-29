# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -Wextra -O3

# Directories
SRC_DIR = src
BUILD_DIR = build
EXE_DIR = exe

# Source files and object files
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS))

# Output executable
TARGET = $(EXE_DIR)/pn_nbody

# Default target
all: $(BUILD_DIR) $(EXE_DIR) $(TARGET)

# Rule to build the target
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ -L/Users/fheinze/Desktop/Cuba-4.2.2 -lcuba -lm

# Rule to build object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Ensure the build directory exists
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Ensure the exe directory exists
$(EXE_DIR):
	mkdir -p $(EXE_DIR)

# Clean up build artifacts
clean:
	rm -rf $(BUILD_DIR) $(EXE_DIR)

# Phony targets
.PHONY: all clean

