CC=gcc
CFLAGS=-O3 -march=native -Wall -idirafter ./src/dsBenchmark # -fno-builtin -g

OBJDIR=obj
SRCDIR=src
OBJS=$(addprefix $(OBJDIR)/,ZigguratO.o aux.o sampler.o permutations.o scheme.o auxzkp.o nizkp.o test_scheme.o)
PROGRAM=tests

ifdef ssl
all: $(OBJS) $(PROGRAM)
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@ -L$(ssl) -lcrypto -lmpfr -lgmp -lflint -lm -fopenmp
$(PROGRAM): $(OBJS)
	@mkdir -p ./benchmarks
	@mkdir -p ./benchmarks/100-512-2
	@mkdir -p ./benchmarks/100-1024-2
	@mkdir -p ./benchmarks/100-1024-4
	$(CC) $(CFLAGS) $(OBJS) -o $(PROGRAM) -L$(ssl) -lcrypto -lmpfr -lgmp -lflint -lm -fopenmp
else
all: $(OBJS) $(PROGRAM)
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@ -lssl -lcrypto -lmpfr -lgmp -lflint -lm -fopenmp
$(PROGRAM): $(OBJS)
	@mkdir -p ./benchmarks
	@mkdir -p ./benchmarks/100-512-2
	@mkdir -p ./benchmarks/100-1024-2
	@mkdir -p ./benchmarks/100-1024-4
	$(CC) $(CFLAGS) $(OBJS) -o $(PROGRAM) -lssl -lcrypto -lmpfr -lgmp -lflint -lm -fopenmp
endif

clean:;
	rm -f $(PROGRAM)
	rm -f $(OBJDIR)/*.o
