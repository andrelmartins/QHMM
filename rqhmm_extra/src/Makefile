SRC:=rqhmm_extra.cpp

all: rqhmm_extra.so

rqhmm_extra.so: $(SRC)
	R CMD SHLIB -o $@ $(SRC)

clean:
	rm -f rqhmm_extra.so $(patsubst %.cpp,%.o,$(SRC))
