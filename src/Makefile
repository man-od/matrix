cc = gcc -Wall -Werror -Wextra -std=c11
GCOVFLAGS = -fprofile-arcs -ftest-coverage
BREAKPOINT = b *.c:12

OS = $(shell uname)
ifeq ($(OS), Linux)
 CHECK_FLAGS = -lcheck -pthread -lrt -lm -lsubunit
else
 CHECK_FLAGS = -lcheck -lm
endif

all: s21_matrix.a

s21_matrix.a:
	$(cc) -c s21_matrix.c
	ar rcs s21_matrix.a *.o
	rm -rf *.o
	ranlib s21_matrix.a

test: clean s21_matrix.a
	$(cc) s21_matrix.a *.c $(CHECK_FLAGS) -o s21_test -L .
	./s21_test


style:
	cp ../materials/linters/.clang-format .
	clang-format -n *.c
	clang-format -n *.h

gcov: clean
	$(cc) -fprofile-arcs -ftest-coverage s21_test.c s21_matrix.c $(CHECK_FLAGS) -o unit_test
	./unit_test
	lcov -t "unit_test" -o unit_test.info --no-external -c -d .
	genhtml -o report unit_test.info
	open ./report/index.html

install_brew:
	cd ~
	curl -fsSL https://rawgit.com/kube/42homebrew/master/install.sh | zsh
	brew install lcov

valgrind: test
	valgrind --trace-children=yes --track-fds=yes --track-origins=yes --leak-check=full --show-leak-kinds=all ./s21_test

leak: test
	CK_FORK=no leaks --atExit -- ./s21_test

clean:
	rm -rf *.o *.gcno *.gcda *.a *.info report .clang-format s21_test *.out unit_test