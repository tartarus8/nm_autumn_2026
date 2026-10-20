package nonlinear

import (
	"fmt"
	"go/ast"
	"go/parser"
	"go/token"
	"math"
	"nm/internal/linalg"
	"strconv"
	"strings"
)

type Expr func([]float64) float64
type Vector func([]float64) []float64
type Jacobian func([]float64) linalg.Matrix

func Parse(source string, n int) (Expr, error) {
	tree, err := parser.ParseExpr(source)
	if err != nil {
		return nil, err
	}
	var compile func(ast.Expr) (Expr, error)
	compile = func(node ast.Expr) (Expr, error) {
		switch node := node.(type) {
		case *ast.BasicLit:
			value, err := strconv.ParseFloat(node.Value, 64)
			if err != nil {
				return nil, err
			}
			return func([]float64) float64 {
				return value
			}, nil
		case *ast.ParenExpr:
			return compile(node.X)
		case *ast.Ident:
			if node.Name == "pi" {
				return func([]float64) float64 {
					return math.Pi
				}, nil
			}
			if node.Name == "e" {
				return func([]float64) float64 {
					return math.E
				}, nil
			}
			name := node.Name
			if name == "x" {
				name = "x1"
			}
			i, err := strconv.Atoi(strings.TrimPrefix(name, "x"))
			if !strings.HasPrefix(name, "x") || err != nil || i < 1 || i > n {
				return nil, fmt.Errorf("unknown variable %q", node.Name)
			}
			return func(x []float64) float64 {
				return x[i-1]
			}, nil
		case *ast.UnaryExpr:
			value, err := compile(node.X)
			if err != nil {
				return nil, err
			}
			if node.Op == token.SUB {
				return func(x []float64) float64 {
					return -value(x)
				}, nil
			}
			if node.Op == token.ADD {
				return value, nil
			}
		case *ast.BinaryExpr:
			left, err := compile(node.X)
			if err != nil {
				return nil, err
			}
			right, err := compile(node.Y)
			if err != nil {
				return nil, err
			}
			switch node.Op {
			case token.ADD:
				return func(x []float64) float64 {
					return left(x) + right(x)
				}, nil
			case token.SUB:
				return func(x []float64) float64 {
					return left(x) - right(x)
				}, nil
			case token.MUL:
				return func(x []float64) float64 {
					return left(x) * right(x)
				}, nil
			case token.QUO:
				return func(x []float64) float64 {
					return left(x) / right(x)
				}, nil
			}
		case *ast.CallExpr:
			name, ok := node.Fun.(*ast.Ident)
			if !ok {
				return nil, fmt.Errorf("only named math functions are allowed")
			}
			if name.Name == "pow" && len(node.Args) == 2 {
				a, err := compile(node.Args[0])
				if err != nil {
					return nil, err
				}
				b, err := compile(node.Args[1])
				if err != nil {
					return nil, err
				}
				return func(x []float64) float64 {
					return math.Pow(a(x), b(x))
				}, nil
			}
			functions := map[string]func(float64) float64{"sqrt": math.Sqrt, "exp": math.Exp, "log": math.Log, "sin": math.Sin, "cos": math.Cos, "tan": math.Tan, "abs": math.Abs}
			fn, ok := functions[name.Name]
			if !ok || len(node.Args) != 1 {
				return nil, fmt.Errorf("unsupported function or argument count: %s", name.Name)
			}
			arg, err := compile(node.Args[0])
			if err != nil {
				return nil, err
			}
			return func(x []float64) float64 {
				return fn(arg(x))
			}, nil
		}
		return nil, fmt.Errorf("unsupported expression in %q; use +, -, *, / and pow for powers", source)
	}
	return compile(tree)
}

func CompileVector(sources []string, n int) (Vector, error) {
	expressions := make([]Expr, len(sources))
	for i, source := range sources {
		e, err := Parse(source, n)
		if err != nil {
			return nil, err
		}
		expressions[i] = e
	}
	return func(x []float64) []float64 {
		values := make([]float64, len(expressions))
		for i, e := range expressions {
			values[i] = e(x)
		}
		return values
	}, nil
}

func CompileJacobian(sources [][]string, n int) (Jacobian, error) {
	rows := make([]Vector, len(sources))
	for i, source := range sources {
		row, err := CompileVector(source, n)
		if err != nil {
			return nil, err
		}
		rows[i] = row
	}
	return func(x []float64) linalg.Matrix {
		matrix := make(linalg.Matrix, len(rows))
		for i, row := range rows {
			matrix[i] = row(x)
		}
		return matrix
	}, nil
}
