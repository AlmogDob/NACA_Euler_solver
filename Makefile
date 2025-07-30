CFLAGS = -Wall -Wextra -std=c99 -lm -lsqlite3 -O3
CCHECKS = -fsanitize=address -g
O_FILES_MAIN = ./build/main.o ./build/mesher.o ./build/solver.o
O_FILES_TEMP = ./build/temp.o ./build/mesher.o ./build/solver.o

# IN_FILE=mesher_input.txt OUT_DIR=./results make main
main: build_and_link_main
	@echo
	./build/main $(IN_FILE) $(OUT_DIR)

	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_MAIN) ./build/main

	@echo
	@echo [INFO] done

build_main: 
	@echo [INFO] building main
	@gcc -c ./src/main.c $(CFLAGS) -o ./build/main.o

build_mesher: ./src/mesher.c
	@echo [INFO] building mesher
	@gcc -c ./src/mesher.c $(CFLAGS) -o ./build/mesher.o

build_solver: ./src/solver.c
	@echo [INFO] building solver
	@gcc -c ./src/solver.c $(CFLAGS) -o ./build/solver.o

link_main: $(O_FILES_MAIN)
	@echo [INFO] linking
	@gcc $(O_FILES_MAIN) $(CFLAGS) -o ./build/main

build_and_link_main: build_mesher build_solver build_main link_main

clean_main:
	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_MAIN) ./build/main


debug_main: debug_build_mesher debug_build_main debug_build_solver link_main
	gdb ./build/main

	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_MAIN)

	@echo
	@echo [INFO] done

debug_build_main: 
	@echo [INFO] building main
	@gcc -c ./src/main.c $(CFLAGS) -ggdb -o ./build/main.o

debug_build_mesher: ./src/mesher.c
	@echo [INFO] building mesher
	@gcc -c ./src/mesher.c $(CFLAGS) -ggdb -o ./build/mesher.o

debug_build_solver: ./src/solver.c
	@echo [INFO] building solver
	@gcc -c ./src/solver.c $(CFLAGS) -ggdb -o ./build/solver.o

profile_main: profile_build_mesher profile_build_solver profile_build_main profile_link_main
	./build/main $(IN_FILE) $(OUT_DIR)
	@echo
	gprof ./build/main gmon.out | /home/almog/.local/bin/gprof2dot | dot -Tpng -Gdpi=200 -o output.png
	# @sleep 0.1
	imview ./output.png
	@echo
	rm ./gmon.out ./output.png 
	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_MAIN)

profile_link_main: $(O_FILES_MAIN)
	@echo [INFO] linking
	@gcc $(O_FILES_MAIN) $(CFLAGS) -p -ggdb -o ./build/main

profile_build_main: 
	@echo [INFO] building main
	@gcc -c ./src/main.c $(CFLAGS) -p -ggdb -o ./build/main.o

profile_build_mesher: ./src/mesher.c
	@echo [INFO] building mesher
	@gcc -c ./src/mesher.c $(CFLAGS) -p -ggdb -o ./build/mesher.o

profile_build_solver: ./src/solver.c
	@echo [INFO] building solver
	@gcc -c ./src/solver.c $(CFLAGS) -p -ggdb -o ./build/solver.o

# valgrind -s --leak-check=full ./build/mesher
# cloc --exclude-lang=JSON,make .

###############################################################

temp: build_mesher build_solver build_temp link_temp
	@echo
	./build/temp $(IN_FILE) $(OUT_DIR)

	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_TEMP) ./build/temp

	@echo
	@echo [INFO] done

build_temp: 
	@echo [INFO] building temp
	@gcc -c ./src/temp.c $(CFLAGS) -o ./build/temp.o

link_temp: $(O_FILES_TEMP)
	@echo [INFO] linking
	@gcc $(O_FILES_TEMP) $(CFLAGS) -o ./build/temp

##########################################################
automat: build_automat
	./build/automat

	@echo
	@echo [INFO] removing automat
	rm -r ./build/automat

	@echo
	@echo [INFO] done

build_automat: ./src/automat.c
	@echo [INFO] building automat
	@gcc ./src/automat.c $(CFLAGS) -o ./build/automat

run_automat:
	@echo
	@./src/automat

clean_automat:
	@echo
	rm ./build/automat
