CFLAGS = -Wall -Wextra -lm -lSDL2 -lSDL2_ttf
CCHECKS = -fsanitize=address

mesher: build_mesher run_mesher clean_mesher  
	@echo [INFO] done

build_mesher: mesher.c
	@echo [INFO] building
	@gcc ./mesher.c $(CFLAGS) -o mesher

run_mesher:
	@echo
	./mesher $(IN_FILE) $(OUT_DIR)

clean_mesher:
	@echo
	rm mesher

debug_build_mesher: mesher.c
	@gcc ./mesher.c $(CFLAGS) -g -o mesher
# valgrind -s --leak-check=full ./mesher
# cloc --exclude-lang=JSON,make .