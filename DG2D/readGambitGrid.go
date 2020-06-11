package DG2D

import (
	"bufio"
	"fmt"
	"io"
	"os"
)

func ReadGambit2d(filename string) {
	var (
		file   *os.File
		err    error
		reader *bufio.Reader
	)
	fmt.Printf("Reading file named: %s\n", filename)
	if file, err = os.Open(filename); err != nil {
		panic(fmt.Errorf("unable to read file %s\n %s", filename, err))
	}
	defer file.Close()
	reader = bufio.NewReader(file)
	var line string
	for {
		line, err = reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				panic(err)
			}
		}
		fmt.Printf(line)
	}
}
