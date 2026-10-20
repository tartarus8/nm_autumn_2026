package main

type Task struct {
	ID       string
	Lower    []float64
	Diagonal []float64
	Upper    []float64
	RHS      []float64
}

var tasks = []Task{
	{
		ID:       "5",
		RHS:      []float64{48, 125, -43, 18, -23},
		Lower:    []float64{-5, -5, -9, 1},
		Diagonal: []float64{8, 22, -11, -15, 7},
		Upper:    []float64{4, 8, 1, 1},
	},
	{
		ID:       "6",
		RHS:      []float64{-58, 161, -114, -90, -55},
		Lower:    []float64{-6, 9, 8, 6},
		Diagonal: []float64{6, 16, -17, 22, -13},
		Upper:    []float64{-5, 9, -3, -8},
	},
}
