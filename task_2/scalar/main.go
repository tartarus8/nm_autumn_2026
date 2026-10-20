package main

import (
	"fmt"
	"math"
	"nm/internal/labio"
	"nm/internal/linalg"
)

func run(extraInfo bool) error {
	t := task
	if err := t.Check(); err != nil {
		return err
	}
	fmt.Println("Lab 2.1, provisional variant 1: 2^x - x^2 - 0.5 = 0")
	fmt.Printf("epsilon = %g; initial = %g; reference = %.15g\n", t.Epsilon, t.Initial, t.Reference)
	if extraInfo {
		fmt.Printf("The graph in initial_guess.png crosses zero near 1.574; choose x0 = %g.\n", t.Initial)
		fmt.Printf("phi(x) = sqrt(2^x - 0.5); phi maps [%g, %g] into itself; |phi'(x)| <= q = %g.\n", t.Lower, t.Upper, t.Q)
		fmt.Println("Newton: x_next = x - f(x)/f'(x).")
	}
	counts := []int{}
	for _, method := range []string{"iteration", "newton"} {
		var history []linalg.Step
		var err error
		metric := "error bound"
		if method == "iteration" {
			history, err = iteration(t)
		} else {
			history, err = newton(t)
			metric = "step size"
		}
		if err != nil {
			return err
		}
		last := history[len(history)-1]
		x := last.X[0]
		if x <= 0 {
			return fmt.Errorf("%s found a non-positive root", method)
		}
		counts = append(counts, last.Iteration)
		fmt.Printf("%s: x = %.15g; iterations = %d\n", method, x, last.Iteration)
		fmt.Printf("%s = %.3e; residual = %.3e; reference error = %.3e\n", metric, last.Estimate, last.Residual, math.Abs(x-t.Reference))
		if extraInfo {
			labio.PrintHistory(method+" (estimate: "+metric+")", history, []float64{t.Reference})
		}
	}
	fmt.Printf("Simple iteration: %d iterations; Newton: %d iterations.\n", counts[0], counts[1])
	return nil
}

func main() { labio.Main(run) }
