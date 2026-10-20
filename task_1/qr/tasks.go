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
		Parameters: labio.Parameters{Epsilon: 1e-06, MaxIterations: 10000},
		ID:         "5",
		Matrix: linalg.Matrix{
			{5, 8, -2},
			{7, -2, -4},
			{5, 8, -1},
		},
	},
	{
		Parameters: labio.Parameters{Epsilon: 1e-06, MaxIterations: 10000},
		ID:         "6",
		Matrix: linalg.Matrix{
			{8, -1, -3},
			{-5, 9, -8},
			{4, -5, 7},
		},
	},
	{
		Parameters: labio.Parameters{Epsilon: 1e-12, MaxIterations: 10000},
		ID:         "7",
		Matrix: linalg.Matrix{
			{0, -1, 0},
			{1, 0, 0},
			{0, 0, 2},
		},
	},
}
