package main

import (
	"nm/internal/linalg"
)

type Task struct {
	ID     int
	Matrix linalg.Matrix
	RHS    []float64
}

var tasks = []Task{
	{
		ID: 5,
		Matrix: linalg.Matrix{
			{3, -8, 1, -7},
			{6, 4, 8, 5},
			{-1, 1, -9, -3},
			{-6, 6, 9, -4},
		},
		RHS: []float64{96, -13, -54, 82},
	},
	{
		ID: 6,
		Matrix: linalg.Matrix{
			{1, 2, -1, -7},
			{8, 0, -9, -3},
			{2, -3, 7, 1},
			{1, -5, -6, 8},
		},
		RHS: []float64{-23, 39, -7, 30},
	},
}
