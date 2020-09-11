package utils

import (
	"image/color"
	"time"

	graphics2D "github.com/notargets/avs/geometry"
)

type ColorName uint8

const (
	White ColorName = iota
	Blue
	Red
	Green
	Black
)

func GetColor(name ColorName) (c color.RGBA) {
	switch name {
	case White:
		c = color.RGBA{
			R: 255,
			G: 255,
			B: 255,
			A: 0,
		}
	case Blue:
		c = color.RGBA{
			R: 50,
			G: 0,
			B: 255,
			A: 0,
		}
	case Red:
		c = color.RGBA{
			R: 255,
			G: 0,
			B: 50,
			A: 0,
		}
	case Green:
		c = color.RGBA{
			R: 25,
			G: 255,
			B: 25,
			A: 0,
		}
	case Black:
		c = color.RGBA{
			R: 0,
			G: 0,
			B: 0,
			A: 0,
		}
	}
	return
}

func SleepForever() {
	var ticks int
	for {
		ticks++
		time.Sleep(time.Second)
		if ticks > 1000 {
			break
		}
	}
}

func ArraysTo2Vector(r1, r2 []float64, scaleO ...float64) (g [][2]float64) {
	var (
		scale float64 = 1
	)
	g = make([][2]float64, len(r1))
	if len(scaleO) > 0 {
		scale = scaleO[0]
	}
	for i := range r1 {
		g[i][0] = r1[i] * scale
		g[i][1] = r2[i] * scale
	}
	return
}

func ArraysToPoints(r1, r2 []float64) (points []graphics2D.Point) {
	points = make([]graphics2D.Point, len(r1))
	for i := range r1 {
		points[i].X[0] = float32(r1[i])
		points[i].X[1] = float32(r2[i])
	}
	return
}
