package main

import (
	"math"
	"nm/internal/labio"
	"nm/internal/linalg"
)

type Task struct {
	labio.Parameters
	Initial, Lower, Upper, Reference []float64
	Q                                float64
}

var task = Task{
	Parameters: labio.Parameters{Epsilon: 1e-10, MaxIterations: 1000},
	Initial:    []float64{3, 0.65},
	Lower:      []float64{2.8, 0.6},
	Upper:      []float64{3, 0.7},
	Q:          0.33,
	Reference:  []float64{2.964631378453523, 0.6255356535051563},
}

func f(x []float64) []float64 {
	u, v := x[0], x[1]
	return []float64{(u*u+4)*v - 8, (u-1)*(u-1) + (v-1)*(v-1) - 4}
}

func jacobian(x []float64) linalg.Matrix {
	u, v := x[0], x[1]
	return linalg.Matrix{
		{2 * u * v, u*u + 4},
		{2 * (u - 1), 2 * (v - 1)},
	}
}

func phi(x []float64) []float64 {
	u, v := x[0], x[1]
	return []float64{1 + math.Sqrt(4-(v-1)*(v-1)), 8 / (u*u + 4)}
}
