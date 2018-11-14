CC = mpicxx

PREFIX = /home/tomswinburne/.local

INC = -I${PREFIX}/include -I../src
LIB = -L${PREFIX}/lib

CPPFLAGS = -std=c++11 -O3 -DBOOST_ALL_DYN_LINK
LDFLAGS = -llammps_mpi -lpthread -lboost_random -lboost_system

SRC_DIR = ../src
SRC_FILES = $(wildcard $(SRC_DIR)/*cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,./%.o,$(SRC_FILES))

EXE = pafi

all:
	${CC} ${CPPFLAGS} ${INC} -c ${SRC_FILES}
	${CC} ${CPPFLAGS} ${LIB} ${OBJ_FILES} ${LDFLAGS} -o ${EXE}

clean:
	rm -rf ${OBJ_FILES} ${EXE}
