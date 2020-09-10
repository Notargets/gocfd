package utils

import (
	"image/color"
	"time"
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
