# Compiler
CC = gcc

# Directories
SRC_DIR   := src
BUILD_DIR := build
EXE_DIR   := exe

# Target
TARGET := $(EXE_DIR)/pn-nbody

# Sources / objects
SRCS := $(wildcard $(SRC_DIR)/*.c)
OBJS := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

# Cuba configuration
CUBA_DIR ?=

# If CUBA_DIR is set, use it. Otherwise rely on system paths.
ifeq ($(strip $(CUBA_DIR)),)
  CUBA_CPPFLAGS :=
  CUBA_LDFLAGS  :=
else
  CUBA_CPPFLAGS := -I$(CUBA_DIR)/include
  CUBA_LDFLAGS  := -L$(CUBA_DIR)
endif

# Libraries
LDLIBS  += -lcuba -lm
LDFLAGS += $(CUBA_LDFLAGS)
CPPFLAGS += $(CUBA_CPPFLAGS)

# Compile flags
CFLAGS ?= -O3
CFLAGS += -Wall -Wextra -MMD -MP

# Default target
all: $(BUILD_DIR) $(EXE_DIR) $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

# Compile
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Directories
$(BUILD_DIR):
	mkdir -p $@

$(EXE_DIR):
	mkdir -p $@

# Include auto-generated header dependencies
-include $(DEPS)

# Convenience targets
debug: CFLAGS := -O0 -g -Wall -Wextra -MMD -MP
debug: all

clean:
	rm -rf $(BUILD_DIR) $(EXE_DIR)

.PHONY: all clean debug
