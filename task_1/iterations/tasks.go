package main

import (
	"nm/internal/labio"
	"nm/internal/linalg"
)

type Task struct {
	labio.Parameters
	ID     string
	Matrix linalg.Matrix
	RHS    []float64
}

var tasks = []Task{
	{
		Parameters: labio.Parameters{Epsilon: 0.01, MaxIterations: 10000},
		ID:         "5",
		Matrix: linalg.Matrix{
			{20, 5, 7, 1},
			{-1, 13, 0, -7},
			{4, -6, 17, 5},
			{-9, 8, 4, -25},
		},
		RHS: []float64{-117, -1, 49, -21},
	},
	{
		Parameters: labio.Parameters{Epsilon: 0.01, MaxIterations: 10000},
		ID:         "6",
		Matrix: linalg.Matrix{
			{23, -6, -5, 9},
			{8, 22, -2, 5},
			{7, -6, 18, -1},
			{3, 5, 5, -19},
		},
		RHS: []float64{232, -82, 202, -57},
	},
}
