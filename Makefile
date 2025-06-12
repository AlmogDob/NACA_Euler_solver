CFLAGS = -Wall -Wextra -lm -std=c99
CCHECKS = -fsanitize=address

# IN_FILE=mesher_input.txt make main
main: build_mesher build_main 
	@echo [INFO] linking
	@gcc ./build/main.o ./build/mesher.o $(CFLAGS) -o ./build/main

	@echo
	./build/main $(IN_FILE)

	@echo
	@echo [INFO] removing build files
	rm -r ./build/main ./build/main.o ./build/mesher.o

	@echo
	@echo [INFO] done

build_main: 
	@echo [INFO] building main
	@gcc -c ./src/main.c $(CFLAGS) -o ./build/main.o

build_mesher: ./src/mesher.c
	@echo [INFO] building mesher
	@gcc -c ./src/mesher.c $(CFLAGS) -o ./build/mesher.o


debug_main: debug_build_mesher debug_build_main
	@echo [INFO] linking
	@gcc ./build/main.o ./build/mesher.o $(CFLAGS) -o ./build/main

	gdb ./build/main

	@echo
	@echo [INFO] removing build files
	rm -r ./build/main ./build/main.o ./build/mesher.o

	@echo
	@echo [INFO] done

debug_build_main: 
	@echo [INFO] building main
	@gcc -c ./src/main.c $(CFLAGS) -ggdb -o ./build/main.o

debug_build_mesher: ./src/mesher.c
	@echo [INFO] building mesher
	@gcc -c ./src/mesher.c $(CFLAGS) -ggdb -o ./build/mesher.o


# valgrind -s --leak-check=full ./build/mesher
# cloc --exclude-lang=JSON,make .