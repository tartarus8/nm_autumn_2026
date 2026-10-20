package main

import (
	"math"
	"nm/internal/labio"
)

type Task struct {
	labio.Parameters
	Initial, Lower, Upper, Q, Reference float64
}

var task = Task{
	Parameters: labio.Parameters{Epsilon: 1e-10, MaxIterations: 1000},
	Initial:    1.6,
	Lower:      1.4,
	Upper:      1.6,
	Q:          0.67,
	Reference:  1.573828923644949,
}

func f(x float64) float64 {
	return math.Pow(2, x) - x*x - 0.5
}

func derivative(x float64) float64 {
	return math.Log(2)*math.Pow(2, x) - 2*x
}

func phi(x float64) float64 {
	return math.Sqrt(math.Pow(2, x) - 0.5)
}
