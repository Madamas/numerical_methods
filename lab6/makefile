TARGET = lab6

.PHONY = all clean

all: $(TARGET)

clean: 
	rm -rvf $(TARGET)
	rm -rvf table.csv

$(TARGET): $(TARGET).cpp
	g++ -o $(TARGET) $(TARGET).cpp