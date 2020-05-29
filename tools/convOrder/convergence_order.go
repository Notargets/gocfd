package main

import (
	"bufio"
	"encoding/csv"
	"flag"
	"fmt"
	"os"
	"strconv"
)

var (
	csvFile string
)

func main() {
	csvFilePtr := flag.String("csvFile", csvFile, "file containing entries of a convergence study")
	flag.Parse()
	csvFile = *csvFilePtr
	if len(csvFile) == 0 {
		flag.Usage()
		os.Exit(1)
	}
	fmt.Printf("Input file: %v\n", csvFile)
	studies := readCSV(csvFile)
	for _, cs := range studies {
		fmt.Printf("Title = %s, Order = %d, CFL = %5.2f\n", cs.title, cs.order, cs.CFL)
		for i := range cs.numPTS {
			fmt.Printf("%d, %v, %v, %v, %v, %v, %v\n",
				cs.numPTS[i], cs.rhoRMS[i], cs.rhouRMS[i], cs.eRMS[i], cs.rhoMAX[i], cs.rhouMAX[i], cs.eMAX[i])
		}
	}
}

type ConvergenceStudy struct {
	title                 string
	order                 int
	numPTS                []int
	CFL                   float64
	rhoRMS, rhouRMS, eRMS []float64
	rhoMAX, rhouMAX, eMAX []float64
}

func NewConvergenceStudy(title string, order int, CFL float64) *ConvergenceStudy {
	return &ConvergenceStudy{
		title: title,
		order: order,
		CFL:   CFL,
	}
}

func (cs *ConvergenceStudy) Add(numPTS int, rhoRMS, rhouRMS, eRMS, rhoMAX, rhouMAX, eMAX float64) {
	cs.numPTS = append(cs.numPTS, numPTS)
	cs.rhoRMS = append(cs.rhoRMS, rhoRMS)
	cs.rhouRMS = append(cs.rhouRMS, rhouRMS)
	cs.eRMS = append(cs.eRMS, eRMS)
	cs.rhoMAX = append(cs.rhoMAX, rhoMAX)
	cs.rhouMAX = append(cs.rhouMAX, rhouMAX)
	cs.eMAX = append(cs.eMAX, eMAX)
}

func readCSV(csvFile string) (studies map[string]*ConvergenceStudy) {
	var (
		records                                      [][]string
		err                                          error
		f                                            *os.File
		ok                                           bool
		cs                                           *ConvergenceStudy
		cfl                                          float64
		rhoRMS, rhouRMS, eRMS, rhoMAX, rhouMAX, eMAX float64
	)
	studies = make(map[string]*ConvergenceStudy)
	if f, err = os.Open(csvFile); err != nil {
		panic(err)
	}
	r := csv.NewReader(bufio.NewReader(f))
	if records, err = r.ReadAll(); err != nil {
		panic(err)
	}
	for i, rec := range records {
		if i == 0 {
			continue
		}
		title, nptstxt, ntxt, cfltxt := rec[0], rec[1], rec[2], rec[3]
		n, _ := strconv.Atoi(ntxt)
		npts, _ := strconv.Atoi(nptstxt)
		_, _ = fmt.Sscanf(cfltxt, "%f", &cfl)
		combTitle := title + ntxt
		if cs, ok = studies[combTitle]; !ok {
			cs = NewConvergenceStudy(title, n, cfl)
			studies[combTitle] = cs
		}
		_, _ = fmt.Sscanf(rec[4], "%f", &rhoRMS)
		_, _ = fmt.Sscanf(rec[5], "%f", &rhouRMS)
		_, _ = fmt.Sscanf(rec[6], "%f", &eRMS)
		_, _ = fmt.Sscanf(rec[7], "%f", &rhoMAX)
		_, _ = fmt.Sscanf(rec[8], "%f", &rhouMAX)
		_, _ = fmt.Sscanf(rec[9], "%f", &eMAX)
		cs.Add(npts, rhoRMS, rhouRMS, eRMS, rhoMAX, rhouMAX, eMAX)
	}
	return
}
