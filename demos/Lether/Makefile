.PHONY: all clean

LAZER_DIR = ../..

LIBS = $(LAZER_DIR)/liblazer.a -lmpfr -lgmp -lm ../../third_party/hexl-development/build/hexl/lib/libhexl.a -lstdc++ #XXX
CFLAGS = -Wall -Wextra -Wshadow -Wundef -O3 -g

all: Lether-demo

Lether-demo: params.h Lether-demo.c
	$(CC) $(CFLAGS) -I$(LAZER_DIR) -o $@ Lether-demo.c $(LIBS)


clean:
	rm -f Lether-demo