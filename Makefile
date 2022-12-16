# Write this at some point today please
FILES = min_polynomial finite_field list

CC = clang
CFLAGS = -Iinclude $(shell sdl2-config --cflags | sed -e "s/include\/SDL2/include/") -Wall -g -fno-omit-frame-pointer -fsanitize=address -Wno-nullability-completeness

FILES_OBJ = $(addprefix out/,$(FILES:=.o))

all: $(FILES_OBJ)