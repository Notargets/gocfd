package readfiles

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"github.com/notargets/gocfd/utils"
)

func getToken(reader *bufio.Reader) (token string) {
	var (
		line string
		err  error
	)
	line = getLineNoComments(reader)
	ind := strings.Index(line, "=")
	if ind < 0 {
		err = fmt.Errorf("badly formed input line [%s], should have an =", line)
		panic(err)
	}
	token = line[ind+1:]
	return
}

func readLabel(reader *bufio.Reader) (label string) {
	var (
		err error
	)
	token := getToken(reader)
	if _, err = fmt.Sscanf(token, "%s", &label); err != nil {
		err = fmt.Errorf("unable to read label from token: [%s]", token)
		panic(err)
	}
	label = strings.Trim(label, " ")
	return
}

func readNumber(reader *bufio.Reader) (num int) {
	var (
		err error
	)
	token := getToken(reader)
	if _, err = fmt.Sscanf(token, "%d", &num); err != nil {
		err = fmt.Errorf("unable to read number from token: [%s]", token)
		panic(err)
	}
	return
}

func getLineNoComments(reader *bufio.Reader) (line string) {
	var ()
	for {
		line = strings.Trim(getLine(reader), " ")
		//fmt.Printf("line = [%s]\n", line)
		ind := strings.Index(line, "%")
		if ind < 0 || ind != 0 {
			return
		}
	}
}

func ReadSU2(filename string, verbose bool) (K int, VX, VY utils.Vector, EToV, BCType utils.Matrix) {
	var (
		file   *os.File
		err    error
		reader *bufio.Reader
	)
	if verbose {
		fmt.Printf("Reading SU2 file named: %s\n", filename)
	}
	if file, err = os.Open(filename); err != nil {
		panic(fmt.Errorf("unable to open file %s\n %s", filename, err))
	}
	defer file.Close()
	reader = bufio.NewReader(file)

	dimensionality := readNumber(reader)
	fmt.Printf("Read file with %d dimensional data...\n", dimensionality)

	return
}
