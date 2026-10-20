package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"nm/internal/linalg"
	"os"
	"path/filepath"
	"runtime"
)

type Task struct {
	ID     int           `json:"id"`
	Matrix linalg.Matrix `json:"matrix"`
	RHS    []float64     `json:"rhs"`
}

func inputPath() string {
	_, file, _, _ := runtime.Caller(0)
	return filepath.Join(filepath.Dir(file), "tests.json")
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

func run(path string, extraInfo bool) error {
	tasks, err := readTasks(path)
	if err != nil {
		return err
	}
	for _, task := range tasks {
		fmt.Printf("Task %d\n\n", task.ID)
		if extraInfo {
			linalg.PrintMatrix("A", task.Matrix)
			fmt.Printf("b = %.10g\n\nLU factorization\n", task.RHS)
		}
		f, err := linalg.FactorLU(task.Matrix, extraInfo)
		if err != nil {
			return fmt.Errorf("task %d: %w", task.ID, err)
		}
		if extraInfo {
			fmt.Println("Solve Ax = b")
		}
		x := linalg.SolveLU(f, task.RHS, extraInfo)
		det := linalg.DeterminantLU(f)
		if extraInfo {
			fmt.Printf("\ndet(A) = (-1)^%d * product(diag(U)) = %.10g\n", f.Swaps, det)
		}
		inverse := linalg.InverseLU(f, extraInfo)

		fmt.Printf("\nx = %.10g\ndet(A) = %.10g\n\n", x, det)
		linalg.PrintMatrix("L", f.L)
		linalg.PrintMatrix("U", f.U)
		linalg.PrintMatrix("P", f.P)
		linalg.PrintMatrix("A^-1", inverse)

		lu := linalg.Multiply(f.L, f.U)
		pa := linalg.Multiply(f.P, task.Matrix)
		product := linalg.Multiply(task.Matrix, inverse)
		ax := linalg.MatVec(task.Matrix, x)
		linalg.PrintMatrix("L * U", lu)
		linalg.PrintMatrix("P * A", pa)
		fmt.Printf("max |LU - PA| = %.3e\n\n", linalg.MaxDifference(lu, pa))
		fmt.Printf("Ax = %.10g\nb  = %.10g\n", ax, task.RHS)
		fmt.Printf("max |Ax - b| = %.3e\n\n", linalg.MaxDifference(linalg.Matrix{ax}, linalg.Matrix{task.RHS}))
		linalg.PrintMatrix("A * A^-1", product)
		fmt.Printf("max |A * A^-1 - I| = %.3e\n\n", linalg.MaxDifference(product, linalg.Identity(len(task.Matrix))))
	}
	return nil
}

func main() {
	extraInfo := flag.Bool("extraInfo", false, "print pivoting, elimination, and substitution steps")
	flag.Parse()
	if err := run(inputPath(), *extraInfo); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}
