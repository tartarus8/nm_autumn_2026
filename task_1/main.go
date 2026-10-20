package main

import (
	"encoding/json"
	"fmt"
	"os"
)

type Task struct {
	ID     int       `json:"id"`
	Matrix Matrix    `json:"matrix"`
	RHS    []float64 `json:"rhs"`
}

func readTasks(path string) ([]Task, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}

	var input struct {
		Tasks []Task `json:"tasks"`
	}

	if err := json.Unmarshal(data, &input); err != nil {
		return nil, err
	}

	for _, task := range input.Tasks {
		if len(task.Matrix) == 0 {
			return nil, fmt.Errorf("task %d: matrix must be nonempty", task.ID)
		}
		for _, row := range task.Matrix {
			if len(row) != len(task.Matrix) {
				return nil, fmt.Errorf("task %d: matrix must be square", task.ID)
			}
		}
		if len(task.RHS) != len(task.Matrix) {
			return nil, fmt.Errorf("task %d: right-hand side length must match the matrix dimension", task.ID)
		}
	}

	return input.Tasks, nil
}

func printMatrix(name string, a Matrix) {
	fmt.Println(name)
	for _, row := range a {
		for _, value := range row {
			fmt.Printf("%12.6g", value)
		}
		fmt.Println()
	}
	fmt.Println()
}

func run(path string, extraInfo bool) error {
	tasks, err := readTasks(path)
	if err != nil {
		return err
	}

	for _, task := range tasks {
		fmt.Printf("Task %d", task.ID)

		if extraInfo {
			printMatrix("A", task.Matrix)
			fmt.Printf("b = %.10g\n\nLU factorization\n", task.RHS)
		}

		f, err := factorLU(task.Matrix, extraInfo)
		if err != nil {
			return fmt.Errorf("task %d: %w", task.ID, err)
		}

		if extraInfo {
			fmt.Println("Solve Ax = b")
		}

		x := solveLU(f, task.RHS, extraInfo)
		det := determinantLU(f)

		if extraInfo {
			fmt.Printf("\ndet(A) = (-1)^%d * product(diag(U)) = %.10g", f.Swaps, det)
		}

		inverse := inverseLU(f, extraInfo)

		fmt.Printf("\nx = %.10g\ndet(A) = %.10g\n", x, det)
		printMatrix("A^-1", inverse)

		if extraInfo {
			lu := multiply(f.L, f.U)
			pa := multiply(f.P, task.Matrix)
			product := multiply(task.Matrix, inverse)
			ax := matVec(task.Matrix, x)

			printMatrix("L * U", lu)
			printMatrix("P * A", pa)
			fmt.Printf("max |LU - PA| = %.3e\n", maxDifference(lu, pa))
			fmt.Printf("Ax = %.10g\nb  = %.10g\n", ax, task.RHS)
			fmt.Printf("max |Ax - b| = %.3e\n", maxDifference(Matrix{ax}, Matrix{task.RHS}))
			printMatrix("A * A^-1", product)
			fmt.Printf("max |A * A^-1 - I| = %.3e\n", maxDifference(product, identity(len(task.Matrix))))
		}
	}

	return nil
}

func main() {
	extraInfo := true
	path := "tests.json"

	if err := run(path, extraInfo); err != nil {
		fmt.Println(err)
	}
}
