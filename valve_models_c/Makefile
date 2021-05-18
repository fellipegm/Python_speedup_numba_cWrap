LINK_TARGET=./bin/valve_sim

CFLAGS+=-lm -Wall

# INCDIR = ./inc
# CFLAGS += -I$(INCDIR)

.PHONY: all

all: CFLAGS+= -O2
all: clean $(LINK_TARGET)

.PHONY: debug
debug: CFLAGS+= -O0 -g
debug: clean $(LINK_TARGET)

# Delete all build output
clean: create
	rm -rf ./bin/* ./obj/*

# Create build output directories if they don't exist
create:
	mkdir -p bin obj

./obj/test_models.o: ./test_models.c
	$(CC) -o $@ -c $^ $(CFLAGS)

./obj/csv_utils.o: ./csv_utils.c
	$(CC) -o $@ -c $^ $(CFLAGS)

./obj/valve_models.o: ./valve_models.c
	$(CC) -o $@ -c $^ $(CFLAGS)

./bin/valve_sim: ./obj/test_models.o ./obj/csv_utils.o ./obj/valve_models.o
	$(CC) -o $@ $^ $(CFLAGS)