.PHONY: clean

TARGET = ./qr
FLAGS = -Wall -Wextra -Wextra -Werror -Wpedantic -g -O3
OBJS = qr.o
DEPS = $(OBJS:.o=.d)

$(TARGET): $(OBJS)
	gcc -o $@ $^ $(LIBS)

-include $(DEPS)

%.o: %.c
	gcc $(FLAGS) -c -o $@ $<
	gcc $(FLAGS) -MM -o $(patsubst %.o, %.d, $@) $<

clean:
	@rm $(OBJS) $(TARGET)
	@find . -name "*.d" | xargs rm
