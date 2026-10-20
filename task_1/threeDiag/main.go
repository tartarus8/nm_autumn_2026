package main

import (
	"fmt"
	"nm/internal/labio"
	"nm/internal/linalg"
)

func run(extraInfo bool) error {
	for _, t := range tasks {
		x, p, q, err := linalg.ThreeDiag(t.Lower, t.Diagonal, t.Upper, t.RHS)
		if err != nil {
			return fmt.Errorf("%s: %w", t.ID, err)
		}
		fmt.Printf("task %s\np = %.10g\nq = %.10g\nx = %.10g\n", t.ID, p, q, x)
		fmt.Println()
	}
	return nil
}
func main() {
	labio.Main(run)
}
