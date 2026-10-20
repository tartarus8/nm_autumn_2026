package main

import (
	"nm/internal/labio"
	"nm/internal/linalg"
)

type Task struct {
	labio.Parameters
	ID     string
	Matrix linalg.Matrix
}

var tasks = []Task{
	{
		Parameters: labio.Parameters{Epsilon: 0.01, MaxIterations: 10000},
		ID:         "5",
		Matrix: linalg.Matrix{
			{0, -7, 7},
			{-7, -9, -5},
			{7, -5, -1},
		},
	},
	{
		Parameters: labio.Parameters{Epsilon: 0.01, MaxIterations: 10000},
		ID:         "6",
		Matrix: linalg.Matrix{
			{9, 2, -7},
			{2, -4, -1},
			{-7, -1, 1},
		},
	},
}
